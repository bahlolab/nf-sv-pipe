
process DELLY_CALL {
    label 'delly'
    label 'C2M4T48'
    tag { sam }

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_idx)
    path(excl)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DELLY.raw.bcf"
    """
    delly call -h ${task.cpus} -q ${params.min_mapq} -g $ref_fa -x $excl -o $out_bcf $bam
    """
}
