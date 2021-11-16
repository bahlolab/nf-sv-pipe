
params.sv_min_size = 500

process filter_duphold {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/filter_duphold", mode: 'symlink'
    tag { fam }

    input:
    tuple val(fam), path(vcf), path(tbi)

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${fam}.filtered.vcf.gz"
    """
    bcftools view $vcf -Ou |
        bcftools filter -Ou -s "DEL_BAD_FC" \\
            -e 'SVTYPE="DEL" & SVLEN<=-$params.sv_min_size & AVG(FMT/DHFFC)>0.7' |
        bcftools filter -Ou -s "DUP_BAD_FC" \\
            -e 'SVTYPE="DUP" & SVLEN>=$params.sv_min_size & AVG(FMT/DHBFC)<1.3' |
        bcftools view -Oz -o $out_vcf &&
        bcftools index -t $out_vcf
    """
}