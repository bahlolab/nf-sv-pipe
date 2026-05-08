
process DELLY_BCF_TO_VCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/delly_bcf_to_vcf", mode: 'symlink'

    input:
    tuple val(sam), path(bcf), path(csi)

    output:
    tuple val(sam), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "${bcf.name.replaceAll(/\\.bcf$/, '')}.vcf.gz"
    """
    bcftools view $bcf -Oz -o $out_vcf
    bcftools index $out_vcf
    """
}
