
process DELLY_CNV_NORM {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }

    input:
    tuple val(sam), path(bcf), path(csi)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.delly_cnv_norm.bcf"
    """
    bcftools view ${bcf} \\
        | delly_cnv_norm.awk \\
        | bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
