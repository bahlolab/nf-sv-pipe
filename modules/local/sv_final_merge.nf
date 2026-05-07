
process SV_FINAL_MERGE {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "${params.outdir}/sv_merged", mode: 'copy'

    input:
        tuple path(jasmine_vcf), path(jasmine_idx),
              path(truvari_vcf), path(truvari_idx)

    output:
        tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${params.id}.sv_merged.vcf.gz"
        """
        bcftools concat -a $jasmine_vcf $truvari_vcf \\
            | bcftools sort \\
            | bgzip > $out_vcf
        tabix $out_vcf
        """
}
