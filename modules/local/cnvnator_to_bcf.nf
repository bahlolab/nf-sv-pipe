
process CNVNATOR_TO_BCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/CNVNATOR/to_vcf", mode: 'symlink'

    input:
    tuple val(sam), path(cnvnator_out)
    path(fai)

    output:
    tuple val(sam), path("${sam}.CNVNATOR.bcf"), path("${sam}.CNVNATOR.bcf.csi")

    script:
    """
    cnvnator2VCF.awk -v prefix=$sam -v sample_name=$sam $cnvnator_out |
        bcftools reheader --fai $fai - |
        bcftools view --threads ${task.cpus} -Ob -o ${sam}.CNVNATOR.bcf
    bcftools index --threads ${task.cpus} ${sam}.CNVNATOR.bcf
    """
}
