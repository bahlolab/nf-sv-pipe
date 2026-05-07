
process CNVNATOR_CALL {
    label 'C2M16T4'
    label 'cnvnator'
    tag { sam }
    maxForks 1
    publishDir "${params.progdir}/CNVNATOR/call", mode: 'symlink'

    input:
    tuple val(sam), path(bam), path(bai), val(chrs)
    path(ref_dir)

    output:
    tuple val(sam), path(out)

    script:
    out = "${sam}_CNVnator.out"
    """
    cnvnator -root ${sam}.root -genome ${params.assembly} -chrom $chrs -tree ${bam}
    cnvnator -root ${sam}.root -genome ${params.assembly} -chrom $chrs -his ${params.bin_size} -d $ref_dir
    cnvnator -root ${sam}.root -stat ${params.bin_size}
    cnvnator -root ${sam}.root -partition ${params.bin_size}
    cnvnator -root ${sam}.root -call ${params.bin_size} > $out
    rm ${sam}.root
    """
}
