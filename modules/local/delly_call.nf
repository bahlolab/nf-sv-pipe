
process DELLY_CALL {
    label 'delly'
    label 'C2M4T24'
    tag { sam }
    publishDir "${params.outdir}/DELLY", mode: 'copy', saveAs: { fn -> publish ? fn : null }

    input:
    tuple val(sam), path(bam), path(bai), val(publish)
    tuple path(ref_fa), path(ref_idx)
    path(excl)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DELLY.bcf"
    """
    delly call -h ${task.cpus} -q ${params.min_mapq} -g $ref_fa -x $excl -o $out_bcf $bam
    """
}
