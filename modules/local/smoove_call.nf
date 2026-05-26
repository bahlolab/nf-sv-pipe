
process SMOOVE_CALL {
    label 'smoove'
    tag { fam }
    cpus   { [2, (bams instanceof List ? bams : [bams]).size() * 2].max() }
    memory { [8.GB, (bams instanceof List ? bams : [bams]).size() * 8.GB].max() }
    time   { 8.h * task.attempt }
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }
    storeDir params.cachedir ? "${params.cachedir}/SMOOVE_CALL" : null


    input:
    tuple val(fam), path(bams), path(bais)
    tuple path(ref_fa), path(ref_idx)
    path exclude

    output:
    tuple val(fam), path(out_final), path("${out_final}.csi")

    script:
    out_tmp   = "${fam}-smoove.genotyped.vcf.gz"
    out_final = "${fam}.SMOOVE.vcf.gz"
    bam_list = (bams instanceof List ? bams : [bams]).join(' ')
    """
    smoove call -x \\
        --name ${fam} \\
        --exclude ${exclude} \\
        --fasta ${ref_fa} \\
        --processes ${task.cpus} \\
        --genotype \\
        ${bam_list}
    
    mv $out_tmp $out_final
    mv ${out_tmp}.csi ${out_final}.csi
    """
}
