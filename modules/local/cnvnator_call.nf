
process CNVNATOR_CALL {
    label 'C2M16T4'
    label 'cnvnator'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/CNVNATOR_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    path(ref_dir)
    val(chrs)

    output:
    tuple val(sam), path(out)

    script:
    out = "${sam}_CNVnator.out"
    chrom_arg = chrs?.trim() ? "-chrom ${chrs}" : ""
    """
    cnvnator -root ${sam}.root -genome ${params.assembly} ${chrom_arg} -tree ${bam}
    cnvnator -root ${sam}.root -genome ${params.assembly} ${chrom_arg} -his ${params.bin_size} -d $ref_dir
    cnvnator -root ${sam}.root -stat ${params.bin_size}
    cnvnator -root ${sam}.root -partition ${params.bin_size}
    cnvnator -root ${sam}.root -call ${params.bin_size} > $out
    rm ${sam}.root
    """
}
