
process VCF_HEADER {
    label 'C1M1T1'
    tag { sam }

    input:
    tuple val(sam), path(vcf_in)
    tuple path(ref), path(fai)

    output:
    tuple path(vcf_out), path("${vcf_out}.csi")

    script:
        vcf_out = "${sam}.quickcnv.vcf.gz"
    """
    bcftools reheader $vcf_in --fai $fai -o tmp.vcf
    bcftools view tmp.vcf --threads ${task.cpus}  -Oz -o $vcf_out
    bcftools index $vcf_out
    """
}
