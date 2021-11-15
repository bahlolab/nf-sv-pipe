
process qdnaseq_call {
    cpus 1
    memory '1 GB'
    time '1 h'
    label 'qdnaseq'
    publishDir "progress/qdnaseq_call", mode: 'symlink'
    tag { "$fam:$sam" }

    input:
        tuple val(fam), val(sam), path(depth_bed), path(ref_fai)

    output:
        tuple val(fam), val(sam), path("${fam}.${sam}.vcf.gz")

    script:
    """
    qdnaseq_call.R $params.assembly $ref_fai $depth_bed ${fam}.${sam}
    """
}

