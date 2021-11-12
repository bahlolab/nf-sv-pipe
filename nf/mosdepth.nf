
params.min_mapq = 37

process mosdepth {
    cpus 2
    memory '4 GB'
    time '1 h'
    tag { "$fam:$sam" }
    publishDir "progress/annotate_id", mode: 'symlink'
    container 'quay.io/biocontainers/mosdepth:0.3.2--h01d7912_0'

    input:
        tuple val(fam), val(sam), path(bam), path(bai), path(bed)

    output:
        tuple val(fam), val(sam), path("${fam}.${sam}.regions.bed.gz")

    script:
    """
    mosdepth ${fam}.${sam} $bam \\
        --fast-mode \\
        --no-per-base \\
        --threads $task.cpus \\
        --by $bed \\
        --mapq $params.min_mapq
    """
}

