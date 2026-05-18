
process MATCHA_MERGE {
    container null
    label 'C4M16T4'
    publishDir "${params.progdir}/matcha_merge", mode: 'copy'

    input:
    path(bcfs)
    path(csis)

    output:
    tuple path("${params.id}.cohort.bcf"), path("${params.id}.cohort.bcf.csi")

    script:
    """
    matcha merge \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o ${params.id}.cohort.bcf \\
        ${bcfs.join(' ')}
    """
}
