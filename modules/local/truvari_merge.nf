
process TRUVARI_MERGE {
    label 'bcftools_truvari'
    label 'C4M16T4'
    tag "${chr ?: 'all'}"

    input:
    tuple path(bcfs), path(csis), val(chr)

    output:
    tuple val(chr), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = chr ? "${params.id}.${chr}.TRUVARI.merge.bcf" : "${params.id}.TRUVARI.merge.bcf"
    def region_arg = chr ? "-r ${chr}" : ""
    def filter_cmd = params.truvari_cohort_filter \
        ? "bcftools view --threads ${task.cpus} -i '${params.truvari_cohort_filter}' -Ob -o ${out_bcf} collapsed.bcf" \
        : "mv collapsed.bcf ${out_bcf}"
    """
    bcftools merge -m none --threads ${task.cpus} ${region_arg} -Oz -o merged.vcf.gz ${bcfs.join(' ')}

    bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \\
        --threads ${task.cpus} -Oz -o sv.vcf.gz merged.vcf.gz
    bcftools index -t --threads ${task.cpus} sv.vcf.gz

    bcftools view -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="INS"' \\
        --threads ${task.cpus} -Oz -o bnd.vcf.gz merged.vcf.gz
    bcftools index -t --threads ${task.cpus} bnd.vcf.gz

    truvari collapse \\
        -i sv.vcf.gz \\
        -o sv_collapsed.vcf \\
        -c removed_sv.vcf.gz \\
        --chain \\
        --sizemax -1 \\
        --refdist ${params.truvari_itvl_refdist} \\
        --pctovl  ${params.truvari_itvl_pctovl} \\
        --pctsize 0.0 \\
        --pctseq 0.0
    bcftools sort -Oz -o sv_collapsed.vcf.gz sv_collapsed.vcf && rm sv_collapsed.vcf

    truvari collapse \\
        -i bnd.vcf.gz \\
        -o bnd_collapsed.vcf\\
        -c removed_bnd.vcf.gz \\
        --chain \\
        --refdist ${params.truvari_bnd_refdist} \\
        --bnddist ${params.truvari_bnddist} \\
        --pctovl 0.0 \\
        --pctsize ${params.truvari_bnd_pctsize} \\
        --pctseq 0.0
    bcftools sort -Oz -o bnd_collapsed.vcf.gz bnd_collapsed.vcf && rm bnd_collapsed.vcf
    
    bcftools index -t --threads ${task.cpus} sv_collapsed.vcf.gz
    bcftools index -t --threads ${task.cpus} bnd_collapsed.vcf.gz
    bcftools concat --allow-overlaps --threads ${task.cpus} \\
        sv_collapsed.vcf.gz bnd_collapsed.vcf.gz -Ob -o collapsed.bcf

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
