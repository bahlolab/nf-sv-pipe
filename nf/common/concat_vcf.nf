
// Remove format to resolve inconsistencies between callers

process concat_vcf {
    label 'C2M2T2'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(vcfs), path(indices)

    output:
    tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${params.id}.combined.vcf.gz"
    """
    bcftools concat $vcfs -Ou --allow-overlaps |
        bcftools annotate --remove FORMAT --threads $task.cpus -Oz -o $out_vcf
    bcftools index -t --threads $task.cpus $out_vcf
    """
}
