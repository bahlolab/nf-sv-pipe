
process DELLY_CNV_GENOTYPE {
    label 'delly'
    label 'C2M8T8'
    tag { sam }

    input:
    tuple val(sam), path(bam), path(bai), path(sites_bcf), path(sites_csi)
    tuple path(ref_fa), path(ref_fai)
    tuple path(map_fa), path(map_gzi), path(map_fai)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.delly_cnv_geno.bcf"
    """
    delly cnv -u -v $sites_bcf -g $ref_fa -m $map_fa -o $out_bcf $bam
    """
}
