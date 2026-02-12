
process MOSDEPTH {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label  'mosdepth'
    tag    "$sample"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.regions.bed.gz")

    script:
    """
    mosdepth \\
        --fast-mode \\
        --no-per-base \\
        --threads $task.cpus \\
        --by $params.bin_size \\
        --mapq 30 \\
        $sample \\
        $bam
    """
}