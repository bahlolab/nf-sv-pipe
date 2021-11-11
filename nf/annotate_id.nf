
process annotate_id {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/annotate_id", mode: 'symlink'

    input:
    tuple val(fam), path(vcf)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = vcf.name.replaceAll('.vcf.gz', '.id.vcf.gz')
    """
    bcftools annotate --set-id '${fam}_%ID' $vcf -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}

