
process DYSGU_CALL {
    label 'dysgu'
    label 'C2M16T8'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/DYSGU_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_idx)

    output:
    tuple val(sam), path(out_vcf)

    script:
    out_vcf = "${sam}.dysgu.vcf.gz"
    """
    dysgu call $ref_fa wd_${sam} $bam \\
        -p ${task.cpus} \\
        --clean \\
        --mq ${params.min_mapq} \\
        --symbolic-sv-size 100 \\
        | bgzip > ${out_vcf}
    """
}
