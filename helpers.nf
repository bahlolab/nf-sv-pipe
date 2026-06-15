
Path path(String filename) {
    file(filename, checkIfExists: true)
}

ArrayList<Map> read_tsv(Path path, List<String> names) {
    path.toFile().readLines().with { lines ->
        lines.each { assert it.split('\t').size() == names.size() }
        lines.collect {
            [names, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

def active_branches() {
    def b = []
    if (params.matcha)  b << 'MATCHA'
    if (params.svdb)    b << 'SVDB'
    if (params.truvari) b << 'TRUVARI'
    b
}

ArrayList<Map> read_caller_manifest(Path tsv) {
    def rows = read_tsv(tsv, ['sample', 'caller', 'path'])
    rows.each { r ->
        if (!(r.caller in params.callers)) {
            error "caller_manifest: caller '${r.caller}' (sample ${r.sample}) is not in params.callers ${params.callers}"
        }
        if (!file(r.path).exists())          error "caller_manifest: path does not exist for (${r.sample}, ${r.caller}): ${r.path}"
        if (!file(r.path + '.csi').exists()) error "caller_manifest: missing .csi index for (${r.sample}, ${r.caller}): ${r.path}.csi"
    }
    def dups = rows.groupBy { [it.sample, it.caller] }.findAll { _k, v -> v.size() > 1 }.keySet()
    if (dups) error "caller_manifest: duplicate (sample, caller) rows: ${dups}"
    rows
}

ArrayList<Map> read_merge_manifest(Path tsv) {
    def rows = read_tsv(tsv, ['sample', 'branch', 'path'])
    def active = active_branches()
    rows.each { r ->
        if (!(r.branch in active)) {
            error "merge_manifest: branch '${r.branch}' (sample ${r.sample}) is not in active branches ${active}"
        }
        if (!file(r.path).exists())          error "merge_manifest: path does not exist for (${r.sample}, ${r.branch}): ${r.path}"
        if (!file(r.path + '.csi').exists()) error "merge_manifest: missing .csi index for (${r.sample}, ${r.branch}): ${r.path}.csi"
    }
    def dups = rows.groupBy { [it.sample, it.branch] }.findAll { _k, v -> v.size() > 1 }.keySet()
    if (dups) error "merge_manifest: duplicate (sample, branch) rows: ${dups}"
    // Full-coverage: every sample must have an entry for every active branch.
    rows*.sample.unique().each { sam ->
        def have = rows.findAll { it.sample == sam }*.branch as Set
        def missing = active.findAll { !(it in have) }
        if (missing) error "merge_manifest: sample '${sam}' is missing entries for active branches ${missing}"
    }
    rows
}

// Group per-caller per-sample VCFs by sample, sorted by params.callers priority.
// Output: [sam, [callers], [bcfs], [csis]] with caller order matching params.callers.
def per_sample_by_caller_priority(vcfs) {
    vcfs.map { caller, sam, bcf, csi -> [groupKey(sam.toString(), params.callers.size()), caller, bcf, csi] }
        .groupTuple(by: 0)
        .map { sam, callers, bcfs, csis ->
            def s = [callers, bcfs, csis].transpose()
                .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                .transpose()
            [sam.target, s[0], s[1], s[2]]
        }
}

// Build the caller cache used in normal mode. Returns a Map with:
//   family_fully_cached: caller -> Set<fam> for which every sample in fam has a cached entry for that caller
//   sam_to_fam:          sample -> family (per params.familial)
def build_caller_cache(List caller_rows) {
    def bams = read_tsv(path(params.bams), ['iid', 'bam'])
    def ped  = read_tsv(path(params.ped),  ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])

    def ped_fam    = ped.collectEntries { [(it.iid): it.fid] }
    def sam_to_fam = bams.collectEntries { row -> [(row.iid): params.familial ? ped_fam[row.iid] : row.iid] }

    def fam_members = [:].withDefault { [] as Set }
    sam_to_fam.each { sam, fam -> fam_members[fam] << sam }

    def family_fully_cached = [:]
    params.callers.each { caller ->
        def cached = caller_rows.findAll { it.caller == caller }*.sample as Set
        family_fully_cached[caller] = fam_members.findAll { _f, members -> members.every { it in cached } }.keySet()
    }

    // Warn about caller manifest entries that won't be used (sample's family is not fully cached for that caller).
    caller_rows.each { r ->
        if (!sam_to_fam.containsKey(r.sample)) {
            log.warn "caller_manifest: sample '${r.sample}' is not in params.bams; entry for ${r.caller} will be ignored"
        } else if (!(sam_to_fam[r.sample] in family_fully_cached[r.caller])) {
            log.warn "caller_manifest: ignoring (${r.sample}, ${r.caller}) because family is not fully cached for ${r.caller}"
        }
    }

    [
        rows:                caller_rows,
        family_fully_cached: family_fully_cached,
        sam_to_fam:          sam_to_fam,
    ]
}

def check_test_fixtures() {
    if (!workflow.profile.contains('test')) return
    [params.bams, params.ped, params.ref_fasta].each { f ->
        if (!file(f).exists()) {
            error """\
            Test fixture missing: ${f}
            Generate fixtures first by running:
                bash test/generate_fixtures.sh
            """.stripIndent()
        }
    }
}

def check_callers() {
    def supported = ['MANTA', 'SMOOVE', 'CNVNATOR', 'DELLY', 'DELLY_CNV', 'DYSGU']
    def unsupported = params.callers - supported
    if (unsupported) {
        error "Unsupported callers in params.callers: ${unsupported}. Supported: ${supported}"
    }
    def dups = params.callers.countBy { it }.findAll { k, n -> n > 1 }.keySet()
    if (dups) {
        error "Duplicate callers in params.callers: ${dups}"
    }
}

def check_apply_filters() {
    def not_called = params.apply_filters - params.callers
    if (not_called) {
        error "params.apply_filters contains callers not in params.callers: ${not_called}"
    }
}

def get_chrs_ch() {
    def eff_chr_prefix = params.chr_prefix != null
        ? params.chr_prefix
        : (params.assembly == 'hg38' ? 'chr' : '')

    def chrs_resolved
    if (params.chrs == null) {
        chrs_resolved = []
    } else if (params.chrs == 'auto') {
        chrs_resolved = ((1..22) + ['X', 'Y']).collect { eff_chr_prefix + it }
    } else if (params.chrs instanceof List) {
        chrs_resolved = params.chrs
    } else {
        error "params.chrs must be null, 'auto', or a List of chromosome names; got: ${params.chrs}"
    }

    Channel.value(chrs_resolved)
}
