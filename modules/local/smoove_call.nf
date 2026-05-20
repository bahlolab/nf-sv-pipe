
process SMOOVE_CALL {
    label 'smoove'
    label 'C2M16T4'
    tag { sam }
    storeDir params.cachedir ? "${params.cachedir}/SMOOVE_CALL" : null

    input:
    tuple val(sam), path(bam), path(bai)
    tuple path(ref_fa), path(ref_fai)
    path(exclude)

    output:
    tuple val(sam), path(out_vcf), path("${out_vcf}.csi")

    script:
    smoove_vcf = "${sam}-smoove.genotyped.vcf.gz"
    out_vcf    = "${sam}.SMOOVE.vcf.gz"
    """
    smoove call $bam \\
        --processes $task.cpus \\
        --outdir . \\
        --exclude $exclude \\
        --name $sam \\
        --fasta $ref_fa \\
        --genotype
        
    mv $smoove_vcf $out_vcf
    mv ${smoove_vcf}.csi ${out_vcf}.csi
    """
}
