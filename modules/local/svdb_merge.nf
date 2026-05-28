
process SVDB_MERGE {
    label 'bcftools_svdb'
    label 'C2M4T24'
    tag "${chr ?: 'all'}"

    input:
    tuple path(bcfs), path(csis), val(chr)

    output:
    tuple val(chr), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = chr ? "${params.id}.${chr}.SVDB.merge.bcf" : "${params.id}.SVDB.merge.bcf"
    def region_arg = chr ? "-r ${chr}" : ""

    // svdb --merge requires VCF.gz + TBI; it also can't handle sample names with dots/slashes,
    // so BCFs are renamed to SAM_<safe> on conversion.
    // svdb adds set/FOUNDBY/svdb_origin/SUPP_VEC/VARID during per-sample collapse (SVDB_COLLAPSE);
    // stripping them here prevents duplicate INFO tags when svdb appends them again at cohort level.
    def svdb_tags = 'INFO/set,INFO/FOUNDBY,INFO/svdb_origin,INFO/SUPP_VEC,INFO/VARID'
    def conv_cmds = bcfs.collect { bcf ->
        def n = "SAM_" + bcf.name.replaceAll(/\..*$/, '').replaceAll(/[^A-Za-z0-9_]/, '_')
        "bcftools annotate ${region_arg} -x ${svdb_tags} -Ou ${bcf.name} | bcftools view --threads ${task.cpus} -Oz -o ${n}.vcf.gz && bcftools index -f -t --threads ${task.cpus} ${n}.vcf.gz"
    }.join('\n    ')
    def safe_vcfs = bcfs.collect { bcf ->
        "SAM_" + bcf.name.replaceAll(/\..*$/, '').replaceAll(/[^A-Za-z0-9_]/, '_') + ".vcf.gz"
    }.join(' ')

    // optional PASS/quality filter on the cohort-merged result
    def filter_out = params.svdb_cohort_filter ? "filtered.bcf" : "merged.bcf"
    def filter_cmd = params.svdb_cohort_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.svdb_cohort_filter}' -Ob -o filtered.bcf merged.bcf"
        : ""

    // svdb cohort merge accumulates N*M per-sample INFO fields; keep only the listed tags
    def info_strip_cmd = params.svdb_info_keep
        ? "bcftools annotate -x '^${params.svdb_info_keep.collect { 'INFO/' + it }.join(',')}' -Ob -o ${out_bcf} ${filter_out}"
        : "mv ${filter_out} ${out_bcf}"
    """
    # convert + strip stale per-sample SVDB tags
    ${conv_cmds}

    # cohort merge
    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_bnd_distance} \\
        --overlap ${params.svdb_overlap} \\
        --vcf ${safe_vcfs} \\
        | bcftools view --threads ${task.cpus} -Ob -o merged.bcf

    # filter then trim INFO bloat
    ${filter_cmd}
    ${info_strip_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
