
process bcftools_merge {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/bcftools_merge", mode: 'symlink'
    tag { fam }

    input:
    tuple val(fam), path(vcfs), path(tbis)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${fam}.merged.vcf.gz"
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