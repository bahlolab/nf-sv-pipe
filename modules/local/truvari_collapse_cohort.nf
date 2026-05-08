
process TRUVARI_COLLAPSE_COHORT {
    label 'truvari'
    label 'C4M16T4'
    publishDir "${params.progdir}/truvari_collapse_cohort", mode: 'symlink'

    input:
    tuple path(vcf), path(tbi)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple path(out_vcf)

    script:
    out_vcf = "${params.id}.cohort_truvari.vcf.gz"
    """
    truvari collapse \\
        --keep first \\
        --reference $ref_fa \\
        --refdist   $params.truvari_cohort_refdist \\
        --pctseq    $params.truvari_cohort_pctseq \\
        --pctsize   $params.truvari_cohort_pctsize \\
        --bnddist   $params.truvari_cohort_bnddist \\
        --input $vcf \\
        --output $out_vcf \\
        --removed-output removed.vcf.gz
    """
}
