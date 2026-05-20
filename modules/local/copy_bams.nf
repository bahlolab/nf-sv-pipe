
process COPY_BAMS {
    label 'C2M2T8'
    container null
    tag { sam }

    input:
    tuple val(fam), val(sam), path(bam), path(bai)

    output:
    tuple val(fam), val(sam), path(bam), path(bai)

    script:
    """
    cp `readlink $bam` tmp && mv tmp $bam
    cp `readlink $bai` tmp && mv tmp $bai
    """
}
