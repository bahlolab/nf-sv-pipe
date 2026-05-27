
process MATCHA_MERGE {
    label 'matcha'
    label 'C4M16T4'
    tag   "${chr ?: 'all'}"

    input:
    tuple path(bcfs), path(csis), val(chr)
    val(chr_set)

    output:
    tuple val(chr), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = chr ? "${params.id}.${chr}.cohort.bcf" : "${params.id}.cohort.bcf"
    def chr_arg     = chr     ? "--chrs ${chr}"        : ""
    def chr_set_arg = chr_set ? "--chr-set ${chr_set}" : ""
    def filter_cmd  = params.matcha_cohort_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.matcha_cohort_filter}' -Ob -o ${out_bcf} merged.bcf"
        : "mv merged.bcf ${out_bcf}"
    """
    matcha merge \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        --missing-to-ref \\
        ${chr_arg} \\
        ${chr_set_arg} \\
        -o merged.bcf \\
        ${bcfs.join(' ')}

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
