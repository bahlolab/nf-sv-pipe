
process SNORM {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    tag    "$sample"
    container null

    input:
    tuple val(sample), path(bed)

    output:
    path "${sample}.bins.rds"    , emit: bins
    path "${sample}.coverage.rds", emit: coverage
    

    script:
    """
    snorm.R $bed $sample
    """
}