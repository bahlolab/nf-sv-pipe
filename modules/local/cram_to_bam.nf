
process CRAM_TO_BAM {
    label 'samtools'
    label 'C2M4T4'
    tag { sam }

    input:
    tuple val(fam), val(sam), path(cram), path(crai)
    tuple path(ref_fa), path(ref_idx)

    output:
    tuple val(fam), val(sam), path("${sam}.bam"), path("${sam}.bam.bai")

    script:
    """
    samtools view -@ ${task.cpus - 1} -b -T $ref_fa -o ${sam}.bam $cram
    samtools index ${sam}.bam
    """
}
