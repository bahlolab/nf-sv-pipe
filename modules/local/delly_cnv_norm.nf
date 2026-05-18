
process DELLY_CNV_NORM {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.progdir}/delly_cnv_norm", mode: 'symlink'

    input:
    tuple val(sam), path(bcf), path(csi)

    output:
    tuple val(sam), path("${sam}.delly_cnv.bcf"), path("${sam}.delly_cnv.bcf.csi")

    script:
    """
    bcftools view ${bcf} | delly_cnv_norm.awk | bcftools view --threads ${task.cpus} -Ob -o ${sam}.delly_cnv.bcf
    bcftools index --threads ${task.cpus} ${sam}.delly_cnv.bcf
    """
}
