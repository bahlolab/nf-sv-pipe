
process MERGE {
    cpus     2
    memory { 16 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    container null

    input:
    path(calls)

    output:
    path('quickcnv.merged.vcf')

    script:
    """
    merge.R quickcnv.merged.vcf $calls
    """
}