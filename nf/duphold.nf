
process duphold {
    cpus 2
    memory '2 GB'
    time '1 h'
    publishDir "progress/duphold", mode: 'symlink'
    tag {"$fam:$sam"}

    input:
    tuple val(fam), val(sam), path(vcf), path(tbi), path(bam), path(bai), path(ref), path(fai)

    output:
    tuple val(fam), val(sam), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = "${fam}.${sam}.duphold.vcf.gz"
    """
    bcftools view $vcf -s $sam -Ou | 
        bcftools view -i 'GT="alt" & (SVTYPE="DEL" | SVTYPE="DUP")' -Ob -o sample_call.bcf &&
        bcftools index sample_call.bcf
    bcftools view $vcf -s $sam -Ou | 
        bcftools view -e 'GT="alt" & (SVTYPE="DEL" | SVTYPE="DUP")' -Ob -o sample_no_call.bcf &&
        bcftools index sample_no_call.bcf
    duphold --threads $task.cpus \\
        --vcf sample_call.bcf \\
        --bam $bam \\
        --fasta $ref \\
        --output duphold.bcf &&
        bcftools index duphold.bcf
    bcftools concat duphold.bcf sample_no_call.bcf -a -Oz -o $out_vcf &&
        bcftools index -t $out_vcf
    """
}

