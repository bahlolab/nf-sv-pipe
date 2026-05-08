
process DELLY_GENOTYPE {
    label 'delly'
    label 'C2M8T2'
    tag { sam }
    publishDir "${params.progdir}/delly_genotype", mode: 'symlink'

    input:
    tuple val(sam), path(bam), path(bai), path(sites_bcf), path(sites_csi)
    tuple path(ref_fa), path(ref_fai)
    path(excl)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.delly_geno.bcf"
    """
    delly call -g $ref_fa -x $excl -v $sites_bcf -o $out_bcf $bam
    """
}
