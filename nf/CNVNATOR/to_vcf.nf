
process to_vcf {
    label 'C1M1T1'
    publishDir "${params.progdir}/CNVNATOR/to_vcf", mode: 'symlink'
    tag { sam }

    input:
        tuple val(sam), path(cnvnator_out), path(ref), path(fai)

    output:
        tuple val(sam), path(vcf), path("${vcf}.csi")

    script:
        vcf = "${sam}.CNVnator.vcf.gz"
        """
        cnvnator2VCF.pl $cnvnator_out $ref -prefix $sam |
            bcftools view -Oz -o tmp.vcf.gz
        echo $sam > sam.txt
        bcftools reheader tmp.vcf.gz --fai $fai --samples sam.txt -o $vcf
        bcftools index $vcf
        rm sam.txt tmp.vcf.gz
        """
}
