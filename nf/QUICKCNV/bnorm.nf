
process BNORM {
    cpus     2
    memory { 16 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    container null
    tag "$shard"

    input:
    tuple val(shard), path(bins), path(snorm)

    output:
    path('*.bnorm.rds')

    script:
    """
    bnorm.R $bins $snorm
    """
}