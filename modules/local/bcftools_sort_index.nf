
process BCFTOOLS_SORT_INDEX {
    label 'bcftools'
    label 'C2M4T4'
    tag { vcf.name }
    publishDir "${params.progdir}/bcftools_sort_index", mode: 'symlink'

    input:
    path(vcf)

    output:
    tuple path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = vcf.name.replaceAll(/\.vcf\.gz$/, '.sorted.vcf.gz')
    """
    bcftools sort -Oz -o sorted.vcf.gz $vcf
    bcftools query -l sorted.vcf.gz | sed 's/^[0-9]*_//' > corrected_samples.txt
    bcftools reheader -s corrected_samples.txt sorted.vcf.gz -o $out_vcf
    bcftools index $out_vcf
    rm sorted.vcf.gz corrected_samples.txt
    """
}
