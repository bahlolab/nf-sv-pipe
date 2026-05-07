
process BCFTOOLS_MERGE_COHORT {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "${params.progdir}/bcftools_merge_cohort", mode: 'symlink'

    input:
    path(vcfs)

    output:
    tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.cohort_truvari_merged.vcf.gz"
    """
    for vcf in $vcfs; do
        bcftools index --tbi \$vcf
    done
    bcftools merge -m none $vcfs -Oz -o $out_vcf
    bcftools index --tbi $out_vcf
    """
}
