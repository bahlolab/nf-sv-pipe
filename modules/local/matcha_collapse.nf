
process MATCHA_COLLAPSE {
    label 'bcftools'
    label 'C2M16T4'
    tag "$sam"

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)
    val(chrs)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.collapsed.bcf"
    def inputs = [callers, bcfs].transpose().collect { caller, bcf -> "${caller}:${bcf}" }.join(' ')
    def chrs_arg = chrs?.trim() ? "--chrs ${chrs}" : ""
    """
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib

    matcha collapse \\
        ${inputs} \\
        ${chrs_arg} \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o ${out_bcf}
    """
}
