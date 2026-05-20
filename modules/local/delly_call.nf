
process DELLY_CALL {
    label 'delly'
    label 'C2M4T16'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/DELLY_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)
    path(excl)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.delly.bcf"
    """
    delly call -h ${task.cpus} -g $ref_fa -x $excl -o $out_bcf $bam
    """
}
