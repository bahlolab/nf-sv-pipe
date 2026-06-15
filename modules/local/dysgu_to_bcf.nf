
process DYSGU_TO_BCF {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.outdir}/DYSGU", mode: 'copy'

    input:
    tuple val(sam), path(vcf)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DYSGU.bcf"
    """
    bcftools view ${vcf} --threads ${task.cpus} -Ob -o ${out_bcf} 
    bcftools index ${out_bcf} --threads ${task.cpus} 
    """
}
