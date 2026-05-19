
process MATCHA_COLLAPSE {
    container null
    label 'C4M16T4'
    tag { sam }

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.collapsed.bcf"
    def inputs = [callers, bcfs].transpose().collect { caller, bcf -> "${caller}:${bcf}" }.join(' ')
    """
    matcha collapse \\
        ${inputs} \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o ${out_bcf}
    """
}
