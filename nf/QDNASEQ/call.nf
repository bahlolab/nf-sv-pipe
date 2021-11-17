
process call {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'qdnaseq'
    publishDir "progress/QDNASEQ/call", mode: 'symlink'
    tag { sam }

    input:
        tuple val(sam), path(depth_bed), path(ref_fai)

    output:
        tuple val(sam), path("${sam}.vcf.gz")

    script:
    """
    qdnaseq_call.R $params.assembly $ref_fai $depth_bed $sam
    """
}

