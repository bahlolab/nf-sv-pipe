
process call {
    cpus 1
    memory '16 GB'
    time '4 h'
    container 'quay.io/biocontainers/cnvnator:0.4.1--h9c7f56d_2'
    publishDir "${params.progdir}/CNVNATOR/call", mode: 'symlink'
    tag { sam }

    input:
        tuple val(sam), path(bam), path(bai), val(chrs), path(ref)

    output:
        tuple val(sam), path(out)

    script:
        out = "${sam}_CNVnator.out"
        """
        cnvnator -root ${sam}.root -genome ${params.assembly} -chrom $chrs -tree ${bam}
        cnvnator -root ${sam}.root -genome ${params.assembly} -chrom $chrs -his ${params.bin_size} -d $ref
        cnvnator -root ${sam}.root -stat ${params.bin_size}
        cnvnator -root ${sam}.root -partition ${params.bin_size}
        cnvnator -root ${sam}.root -call ${params.bin_size} > $out
        rm ${sam}.root
        """
}

/*
    High memory requirement seems to be limited to the -tree module
    This is also the most time consuming so
 */
