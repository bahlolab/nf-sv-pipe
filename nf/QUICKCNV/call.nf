
process CALL {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    tag    "$sample"
    container null
    publishDir "${params.outdir}/QUICKCNV", mode: 'copy'


    input:
    tuple val(sample), path(bnorm)

    output:
    tuple val(sample), path("${sample}.calls.rds"), emit: calls
    tuple val(sample), path("${sample}.bpt.txt")  , emit: reg

    script:
    """
    call.R $sample $bnorm 
    """
}