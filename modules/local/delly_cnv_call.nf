
process DELLY_CNV_CALL {
    label 'delly'
    label 'C2M16T8'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/DELLY_CNV_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_idx)
    tuple path(map_fa), path(map_gzi), path(map_fai)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.delly_cnv.bcf"
    """
    delly cnv -q ${params.min_mapq} -g $ref_fa -m $map_fa -o $out_bcf $bam
    """
}
