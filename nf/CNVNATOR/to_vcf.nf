process to_vcf {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir 'progress/CNVNATOR/to_vcf', mode: 'symlink'

    input:
        tuple val(sam), path(cnvnator_out), path(ref), path(fai)

    output:
        tuple val(sam), path(vcf), path("${vcf}.csi")

    script:
        vcf = "${sam}.CNVnator.vcf.gz"
        """
        cnvnator2VCF.pl $cnvnator_out $ref -prefix $sam |
            bcftools view -Oz -o tmp.vcf.gz
        bcftools reheader tmp.vcf.gz --fai $fai -o $vcf
        bcftools index $vcf
        """
}
