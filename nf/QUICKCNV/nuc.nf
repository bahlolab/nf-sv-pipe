
process NUC {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label  'bedtools'

    input:
    tuple path(ref), path(fai)

    output:
    path('nuc.bed.gz')

    script:
    """
    (
        echo -e "chrom\tstart\tend\tnA\tnC\tnG\tnT\tnN"
        bedtools makewindows -g $fai -w $params.bin_size \\
            | bedtools nuc -fi $ref -bed - \\
            | tail -n+2 \\
            | cut -f1-3,6-10
    ) | gzip > nuc.bed.gz
    """
}