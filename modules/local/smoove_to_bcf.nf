
process SMOOVE_TO_BCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/SMOOVE/to_bcf", mode: 'symlink'

    input:
    tuple val(sam), path(vcf), path(tbi)

    output:
    tuple val(sam), path("${sam}.SMOOVE.bcf"), path("${sam}.SMOOVE.bcf.csi")

    script:
    """
    bcftools view --threads ${task.cpus} -Ob -o ${sam}.SMOOVE.bcf ${vcf}
    bcftools index --threads ${task.cpus} ${sam}.SMOOVE.bcf
    """
}
