
process PASS_FILTER {
    label 'bcftools'
    label 'C2M2T2'
    tag { "${caller}.${sam}" }

    input:
    tuple val(caller), val(sam), path(bcf), path(csi)

    output:
    tuple val(caller), val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.${caller}.pass.bcf"
    """
    bcftools view -f PASS,. ${bcf} --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index ${out_bcf} --threads ${task.cpus}
    """
}
