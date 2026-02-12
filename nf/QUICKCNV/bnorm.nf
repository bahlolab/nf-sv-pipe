
process BNORM {
    cpus     2
    memory { 16 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    container null

    input:
    path(bins)
    path(coverage)

    output:
    path('*.bnorm.rds')

    script:
    """
    bnorm.R $bins $coverage
    """
}