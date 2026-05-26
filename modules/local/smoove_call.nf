
process SMOOVE_CALL {
    label 'smoove'
    tag { fam }
    cpus   { [2, (bams instanceof List ? bams : [bams]).size() * 2].max() }
    memory { [8.GB, (bams instanceof List ? bams : [bams]).size() * 8.GB].max() }
    time   { 8.h * task.attempt }
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    input:
    tuple val(fam), path(bams), path(bais)
    tuple path(ref_fa), path(ref_idx)
    path exclude

    output:
    tuple val(fam), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf  = "${fam}-smoove.genotyped.vcf.gz"
    bam_list = (bams instanceof List ? bams : [bams]).join(' ')
    """
    smoove call -x \\
        --name ${fam} \\
        --exclude ${exclude} \\
        --fasta ${ref_fa} \\
        --processes ${task.cpus} \\
        --genotype \\
        ${bam_list}
    """
}
