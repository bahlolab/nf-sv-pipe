
params.min_mapq = 37
params.exclude_flag = 1540

process mosdepth {
    cpus 2
    memory '4 GB'
    time '1 h'
    tag { sam }
    publishDir "progress/mosdepth", mode: 'symlink'
    container 'quay.io/biocontainers/mosdepth:0.3.2--h01d7912_0'

    input:
        tuple val(sam), path(bam), path(bai), path(bed)

    output:
        tuple val(sam), path("${sam}.regions.bed.gz")

    script:
    """
    mosdepth ${sam} $bam \\
        --fast-mode \\
        --no-per-base \\
        --threads $task.cpus \\
        --by $bed \\
        --mapq $params.min_mapq \\
        --flag $params.exclude_flag
    """
}

