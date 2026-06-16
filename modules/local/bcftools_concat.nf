
process BCFTOOLS_CONCAT {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "$params.outdir", mode: 'copy', saveAs: { fn -> publish ? fn : null }

    input:
    path(bcfs)
    path(csis)
    val tag
    val publish

    output:
    tuple path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${params.id}.${tag}.${publish ? 'cohort' : 'merge'}.bcf"
    """
    bcftools concat --threads ${task.cpus} -Ob -o ${out_bcf} ${bcfs.join(' ')}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
