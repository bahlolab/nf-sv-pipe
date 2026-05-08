
process DELLY_BCF_TO_VCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/delly_bcf_to_vcf", mode: 'symlink'

    input:
    tuple val(sam), path(bcf), path(csi)
    val(is_cnv)

    output:
    tuple val(sam), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = is_cnv ? "${sam}.DELLY_CNV.vcf.gz" : "${sam}.DELLY.vcf.gz"
    """
    ${is_cnv
        ? "bcftools view $bcf | delly_cnv_norm.awk | bcftools view -Oz -o $out_vcf"
        : "bcftools view $bcf -Oz -o $out_vcf"}
    bcftools index $out_vcf
    """
}
