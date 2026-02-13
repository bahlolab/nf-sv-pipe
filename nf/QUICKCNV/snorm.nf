
process SNORM {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    tag    "$sample"
    container null

    input:
    tuple val(sample), path(bed)

    output:
    path "${sample}.shard_*.bins.rds" , emit: bins
    path "${sample}.shard_*.snorm.rds", emit: snorm

    script:
    """
    snorm.R $bed $sample $params.n_phases $params.n_shards
    """
}