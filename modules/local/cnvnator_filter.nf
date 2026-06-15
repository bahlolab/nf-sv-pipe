
process CNVNATOR_FILTER {
    label 'bcftools'
    label 'C2M2T2'
    tag { sam }
    publishDir "${params.outdir}/CNVNATOR", mode: 'copy'

    input:
    tuple val(sam), path(cnvnator_out)
    path(fai)
    path(excl_tsv)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.CNVNATOR.bcf"
    """
    cnvnator2VCF.awk -v prefix=$sam -v sample_name=$sam -v exclude_tsv=$excl_tsv -v min_overlap=${params.cnvnator_exclude_overlap} $cnvnator_out |
        bcftools reheader --fai $fai - |
        bcftools view --threads ${task.cpus} -Ob -o ${out_bcf}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
