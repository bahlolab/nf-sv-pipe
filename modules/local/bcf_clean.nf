
process BCF_CLEAN {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple path(bcf), path(csi)
    val tag
    val info_keep   // List<String>; empty = keep all INFO

    output:
    tuple path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${params.id}.${tag}.cohort.bcf"
    def info_strip = info_keep
        ? "| bcftools annotate -x '^${info_keep.collect { 'INFO/' + it }.join(',')}' -Ou"
        : ""
    """
    bcftools annotate -x '^FORMAT/GT' -Ou ${bcf} \\
        ${info_strip} \\
        | bcftools +setGT -Ou -- -t . -n 0 \\
        | bcftools +fill-tags -Ob -o ${out_bcf} -- -t AN,AC,AF
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
