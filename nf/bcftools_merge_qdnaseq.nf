
process bcftools_merge_qdnaseq {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/bcftools_merge_qdnaseq", mode: 'symlink'
    tag { fam }

    input:
        path(vcfs)

    output:
    tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = params.id + ".qdnaseq.vcf.gz"
    if (vcfs instanceof Path)
        """
        cp $vcfs $out_vcf
        cp $tbis ${out_vcf}.tbi
        """
    else
        """
        bcftools merge -m id $vcfs -Oz -o $out_vcf
        bcftools index -t $out_vcf
        """
}