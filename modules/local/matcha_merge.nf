
process MATCHA_MERGE {
    label 'bcftools'
    label 'C4M16T4'

    input:
    path(bcfs)
    path(csis)

    output:
    tuple path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${params.id}.cohort.bcf"
    """
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib

    matcha merge \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        -o ${out_bcf} \\
        ${bcfs.join(' ')}
    """
}
