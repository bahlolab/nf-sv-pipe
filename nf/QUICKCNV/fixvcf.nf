
process FIXVCF {
    label 'C1M1T1'

    input:
    path(vcf_in)
    tuple path(ref), path(fai)

    output:
    tuple path(vcf_out), path("${vcf_out}.tbi")

    script:
        vcf_out = "quickcnv.vcf.gz"
    """
    bcftools reheader $vcf_in --fai $fai -o tmp.vcf
    bcftools view tmp.vcf --threads ${task.cpus} -Ou \\
        | bcftools +fill-tags -Oz -o $vcf_out -- -t AC,AF,AN
    bcftools index -t $vcf_out
    """
}
