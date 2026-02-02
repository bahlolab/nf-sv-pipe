
process copy_bams {
    label 'C2M2T8'
    publishDir "${params.progdir}/copy_bams", mode: 'symlink'
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
