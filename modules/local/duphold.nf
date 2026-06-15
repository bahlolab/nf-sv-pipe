
process DUPHOLD {
    label 'smoove'
    label 'C2M8T8'
    tag "$sam"
    publishDir { "${params.outdir}/${branch}" }, mode: 'copy'

    input:
    tuple val(branch), val(sam), path(bcf), path(csi), path(bam), path(bai)
    tuple path(ref_fa), path(ref_idx)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.${branch}.duphold.bcf"
    def large_expr = '(INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && (INFO/SVLEN >= ' + params.duphold_min_size + ' || INFO/SVLEN <= -' + params.duphold_min_size + ')'
    def del_cap = params.duphold_max_dels ?
        """
        del_thr=\$(bcftools query -i 'INFO/SVTYPE=="DEL"' -f '[%DHFFC\\n]' large_dh.bcf | sort -g | sed -n '${params.duphold_max_dels}p')
        if [ -n "\$del_thr" ]; then
            del_thr=\$(awk -v a="\$del_thr" -v b="${params.duphold_del_dhffc}" 'BEGIN{print (a<b)?a:b}')
        else
            del_thr=${params.duphold_del_dhffc}
        fi
        """ : "del_thr=${params.duphold_del_dhffc}"
    def dup_cap = params.duphold_max_dups ?
        """
        dup_thr=\$(bcftools query -i 'INFO/SVTYPE=="DUP"' -f '[%DHBFC\\n]' large_dh.bcf | sort -gr | sed -n '${params.duphold_max_dups}p')
        if [ -n "\$dup_thr" ]; then
            dup_thr=\$(awk -v a="\$dup_thr" -v b="${params.duphold_dup_dhbfc}" 'BEGIN{print (a>b)?a:b}')
        else
            dup_thr=${params.duphold_dup_dhbfc}
        fi
        """ : "dup_thr=${params.duphold_dup_dhbfc}"
    """
    bcftools view --threads ${task.cpus} -i '${large_expr}' ${bcf} -Ob -o large.bcf
    bcftools index large.bcf
    bcftools view --threads ${task.cpus} -e '${large_expr}' ${bcf} -Ob -o other.bcf
    bcftools index --threads ${task.cpus} other.bcf

    duphold --threads ${task.cpus} -v large.bcf -b ${bam} -f ${ref_fa} -o large_dh.bcf

    ${del_cap}
    ${dup_cap}
    echo "DUPHOLD ${sam}: del_thr=\$del_thr (default ${params.duphold_del_dhffc}, max_dels ${params.duphold_max_dels ?: 'unset'}); dup_thr=\$dup_thr (default ${params.duphold_dup_dhbfc}, max_dups ${params.duphold_max_dups ?: 'unset'})"

    dh_excl='(INFO/SVTYPE="DEL" && FMT/DHFFC[0] > '\$del_thr') || (INFO/SVTYPE="DUP" && FMT/DHBFC[0] < '\$dup_thr')'
    bcftools view --threads ${task.cpus} -e "\$dh_excl" large_dh.bcf -Ob -o large_pass.bcf
    bcftools index --threads ${task.cpus} large_pass.bcf

    bcftools concat --allow-overlaps --threads ${task.cpus} large_pass.bcf other.bcf -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
