
process DYSGU_CALL {
    label 'dysgu'
    label 'C2M8T8'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/DYSGU_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path(out_vcf)

    script:
    out_vcf = "${sam}.dysgu.vcf.gz"
    """
    dysgu call -p${task.cpus} $ref_fa wd_${sam} $bam \\
        | bgzip > ${out_vcf}
    """
}
