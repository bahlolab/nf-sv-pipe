
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
