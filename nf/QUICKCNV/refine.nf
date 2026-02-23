
process REFINE {
    cpus     2
    memory { 8 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    tag    "$sample"
    container null
    publishDir "${params.outdir}/QUICKCNV", mode: 'copy'


    input:
    tuple val(sample), path(calls), path(bpt_depth)

    output:
    path("${sample}.refined_calls.rds"), emit: calls

    script:
    """
    refine.R $sample $calls $bpt_depth ${params.bin_size*params.n_phases} $params.bin_size
    """
}