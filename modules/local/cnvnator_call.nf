
process CNVNATOR_CALL {
    label 'C2M16T8'
    label 'cnvnator'
    tag { sam }

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
    cnvnator -root ${sam}.root -genome ${params.assembly} ${chrom_arg} -his ${params.cnvnator_bin_size} -d $ref_dir
    cnvnator -root ${sam}.root -stat ${params.cnvnator_bin_size}
    cnvnator -root ${sam}.root -partition ${params.cnvnator_bin_size}
    cnvnator -root ${sam}.root -call ${params.cnvnator_bin_size} > $out
    rm ${sam}.root
    """
}
