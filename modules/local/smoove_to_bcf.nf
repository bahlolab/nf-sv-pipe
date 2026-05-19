
process SMOOVE_TO_BCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }

    input:
    tuple val(sam), path(vcf), path(tbi)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.SMOOVE.bcf"
    """
    bcftools view --threads ${task.cpus} -Ob -o ${out_bcf} ${vcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
