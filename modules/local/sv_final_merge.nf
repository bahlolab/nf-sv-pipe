
process SV_FINAL_MERGE {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "${params.outdir}/sv_merged", mode: 'copy'

    input:
    tuple path(jasmine_vcf), path(jasmine_idx), path(truvari_vcf)

    output:
    tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.sv_merged.vcf.gz"
    """
    bcftools view $truvari_vcf --threads ${task.cpus} --write-index=tbi -Oz -o tmp.truvari.vcf.gz
    bcftools concat $jasmine_vcf tmp.truvari.vcf.gz -a -Oz --write-index=tbi -o $out_vcf
    """
}
