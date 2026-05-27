
process MATCHA_COLLAPSE {
    label 'matcha'
    label 'C2M4T4'
    tag "$sam"
    publishDir "$params.outdir/collapse"

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)
    val(chrs)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.collapsed.bcf"
    def inputs = [callers, bcfs].transpose().collect { caller, bcf -> "${caller}:${bcf}" }.join(' ')
    def chrs_arg = chrs?.trim() ? "--chrs ${chrs}" : ""
    def filter_cmd = params.matcha_sample_filter \
        ? "bcftools view --threads ${task.cpus} -i '${params.matcha_sample_filter}' -Ob -o ${out_bcf} unfiltered.bcf" \
        : "mv unfiltered.bcf ${out_bcf}"
    """
    matcha collapse \\
        ${inputs} \\
        ${chrs_arg} \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o unfiltered.bcf

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
