
process vcf_merge {
    cpus 2
    memory '2 GB'
    time '1 h'
    publishDir "${params.progdir}/MANTA/vcf_merge", mode: 'symlink'

    input:
        path(list)

    output:
        tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${params.id}.duphold.merged.vcf.gz"
        """
        bcftools merge -m id --file-list $list --threads $task.cpus -Oz -o $out_vcf
        bcftools index -t $out_vcf
        """
}