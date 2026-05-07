
process TRUVARI_COLLAPSE_SAMPLE {
    label 'truvari'
    label 'C2M4T4'
    tag { sam }
    publishDir "${params.progdir}/truvari_collapse_sample", mode: 'symlink'

    input:
    tuple val(sam), path(vcf), path(tbi)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path(out_vcf)

    script:
    out_vcf = "${sam}.consensus.vcf.gz"
    """
    truvari collapse \\
        --intra --chain --keep maxqual \\
        --reference $ref_fa \\
        --refdist   $params.truvari_intra_refdist \\
        --pctseq    $params.truvari_intra_pctseq \\
        --pctsize   $params.truvari_intra_pctsize \\
        --bnddist   $params.truvari_intra_bnddist \\
        --input $vcf \\
        --output $out_vcf
    """
}
