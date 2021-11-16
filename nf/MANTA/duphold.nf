
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
        bcftools view -i 'SVTYPE="DEL" | SVTYPE="DUP"' -Ob -o deldup.bcf &&
        bcftools index deldup.bcf
    bcftools view $vcf -s $sam -Ou | 
        bcftools view -e 'SVTYPE="DEL" | SVTYPE="DUP"' -Ob -o other.bcf &&
        bcftools index other.bcf
    duphold --threads $task.cpus \\
        --vcf deldup.bcf \\
        --bam $bam \\
        --fasta $ref \\
        --output duphold.bcf &&
        bcftools index duphold.bcf
    bcftools concat duphold.bcf other.bcf -a -Oz -o $out_vcf &&
        bcftools index -t $out_vcf
    """
}

