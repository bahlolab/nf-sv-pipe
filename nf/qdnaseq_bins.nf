
process qdnaseq_bins {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'qdnaseq'
    publishDir "progress/qdnaseq_bins", mode: 'symlink'

    output:
    path('bins.bed.gz')

    script:
    """
    qdnaseq_bins.R $params.assembly bins.bed.gz
    """
}

