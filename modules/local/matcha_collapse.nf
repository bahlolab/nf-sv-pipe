
process MATCHA_COLLAPSE {
    container null
    label 'C4M16T4'
    tag { sam }
    publishDir "${params.progdir}/matcha_collapse", mode: 'copy'

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)

    output:
    tuple val(sam), path("${sam}.collapsed.bcf"), path("${sam}.collapsed.bcf.csi")

    script:
    def inputs = [callers, bcfs].transpose().collect { caller, bcf -> "${caller}:${bcf}" }.join(' ')
    """
    matcha collapse \\
        ${inputs} \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o ${sam}.collapsed.bcf
    """
}
