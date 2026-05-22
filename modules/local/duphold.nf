
process DUPHOLD {
    label 'smoove'
    label 'C4M8T8'
    tag "$sam"

    input:
    tuple val(sam), path(bcf), path(csi), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.duphold.bcf"
    def large_expr = '(INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP") && (INFO/SVLEN >= ' + params.duphold_min_size + ' || INFO/SVLEN <= -' + params.duphold_min_size + ')'
    def dh_excl   = '(INFO/SVTYPE="DEL" && FMT/DHFFC[0] > ' + params.duphold_del_dhffc + ') || (INFO/SVTYPE="DUP" && FMT/DHBFC[0] < ' + params.duphold_dup_dhbfc + ')'
    """
    bcftools view --threads ${task.cpus} -i '${large_expr}' ${bcf} -Ob -o large.bcf
    bcftools index large.bcf
    bcftools view --threads ${task.cpus} -e '${large_expr}' ${bcf} -Ob -o other.bcf
    bcftools index --threads ${task.cpus} other.bcf

    duphold --threads ${task.cpus} -v large.bcf -b ${bam} -f ${ref_fa} -o large_dh.bcf
    bcftools index --threads ${task.cpus} large_dh.bcf

    bcftools view --threads ${task.cpus} -e '${dh_excl}' large_dh.bcf -Ob -o large_pass.bcf
    bcftools index --threads ${task.cpus} large_pass.bcf

    bcftools concat --allow-overlaps --threads ${task.cpus} large_pass.bcf other.bcf -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
