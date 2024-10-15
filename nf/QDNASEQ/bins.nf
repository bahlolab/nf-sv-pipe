
process bins {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'QDNASEQ'
    publishDir "${params.progdir}/QDNASEQ/bins", mode: 'symlink'

    output:
    path('bins.bed.gz')

    script:
    """
    qdnaseq_bins.R $params.assembly bins.bed.gz
    """
}

