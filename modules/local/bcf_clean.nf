
process BCF_CLEAN {
    label 'bcftools'
    label 'C2M4T4'
    publishDir "$params.outdir", mode: 'copy'

    input:
    tuple path(bcf), path(csi)
    val tag

    output:
    tuple path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${params.id}.${tag}.cohort.bcf"
    """
    bcftools annotate -x '^FORMAT/GT' -Ou ${bcf} \\
        | bcftools +setGT -Ou -- -t . -n 0 \\
        | bcftools +fill-tags -Ob -o ${out_bcf} -- -t AN,AC,AF
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
