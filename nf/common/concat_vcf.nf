

process concat_vcf {
    cpus 2
    memory '1 GB'
    time '1 h'
    publishDir "output", mode: 'copy'

    input:
    tuple path(vcfs), path(indices)

    output:
    tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.merged.vcf.gz"
    """
    bcftools concat $vcfs --threads 2 -Oz -o $out_vcf --allow-overlaps
    bcftools index -t --threads 2 $out_vcf
    """
}
