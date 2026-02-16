
process BPT_DEPTH {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label  'samtools'
    tag    "$sample"

    input:
    tuple val(sample), path(bam), path(bai), path(regions)

    output:
    tuple val(sample), path("${sample}.depth.txt.gz")

    script:
    """
    (
        while read -r REGION; do
            samtools depth $bam \\
            -r \$REGION \\
            -Q 30 \\
            -g SECONDARY \\
            -a \\
            -J
        done < $regions 
    ) | gzip > "${sample}.depth.txt.gz"
    """
}