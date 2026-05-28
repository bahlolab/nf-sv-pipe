
process BCFTOOLS_CONCAT {
    label 'bcftools'
    label 'C2M4T4'

    input:
    path(bcfs)
    path(csis)
    val tag

    output:
    tuple path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${params.id}.${tag}.merge.bcf"
    """
    bcftools concat --threads ${task.cpus} -Ob -o ${out_bcf} ${bcfs.join(' ')}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
