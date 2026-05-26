
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
    def max_calls = params.delly_cnv_max_calls
    def qual_filter = max_calls ?
        """
        min_qual=\$(bcftools query -f '%QUAL\\n' ${bcf} | sort -rn | sed -n '${max_calls}p')
        [ -n "\$min_qual" ] && view_args="-i QUAL>=\$min_qual" || view_args=''
        """ : "view_args=''"
    """
    ${qual_filter}
    bcftools view \$view_args ${bcf} \\
        | delly_cnv_norm.awk \\
        | bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
