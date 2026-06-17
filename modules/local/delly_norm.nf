
process DELLY_NORM {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.outdir}/DELLY", mode: 'copy'

    input:
    tuple val(sam), path(bcf), path(csi)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DELLY.bcf"
    """
    bcftools view ${bcf} | svlen_fix.awk | bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
