
process DELLY_CNV_NORM {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.outdir}/DELLY_CNV", mode: 'copy'

    input:
    tuple val(sam), path(bcf), path(csi)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.DELLY_CNV.bcf"
    def del_cap = params.delly_cnv_max_dels ?
        """
        del_thr=\$(bcftools query -i 'INFO/SVTYPE=="DEL"' -f '%QUAL\\n' normed.bcf | sort -rn | sed -n '${params.delly_cnv_max_dels}p')
        [ -n "\$del_thr" ] && del_excl='INFO/SVTYPE=="DEL" && QUAL<'\$del_thr || del_excl='0'
        """ : "del_excl='0'"
    def dup_cap = params.delly_cnv_max_dups ?
        """
        dup_thr=\$(bcftools query -i 'INFO/SVTYPE=="DUP"' -f '%QUAL\\n' normed.bcf | sort -rn | sed -n '${params.delly_cnv_max_dups}p')
        [ -n "\$dup_thr" ] && dup_excl='INFO/SVTYPE=="DUP" && QUAL<'\$dup_thr || dup_excl='0'
        """ : "dup_excl='0'"
    """
    bcftools view ${bcf} \\
        | delly_cnv_norm.awk \\
        | svlen_fix.awk \\
        | bcftools view --threads ${task.cpus} -Ob -o normed.bcf
    bcftools index --threads ${task.cpus} normed.bcf

    ${del_cap}
    ${dup_cap}
    echo "DELLY_CNV_NORM ${sam}: del_min_qual=\${del_thr:-unset} (max_dels ${params.delly_cnv_max_dels ?: 'unset'}); dup_min_qual=\${dup_thr:-unset} (max_dups ${params.delly_cnv_max_dups ?: 'unset'})"

    bcftools view --threads ${task.cpus} -e "(\$del_excl) || (\$dup_excl)" normed.bcf -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
