
process publish_vcf {
    cpus 2
    memory '1 GB'
    time '1 h'
    publishDir "output", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path(out_vcf), path(out_tbi)

    script:
    out_vcf = "${params.id}.${params.caller}.vcf.gz"
    out_tbi = out_vcf + '.tbi'
    """
    bcftools query -l $vcf | sort > samples.txt
    bcftools view $vcf -Ou -S samples.txt |
        bcftools annotate --threads 2 -Oz -o $out_vcf\\
            --set-id '${params.caller}_%SVTYPE\\_%CHROM\\_%POS\\_%END'
    bcftools index -t $out_vcf --threads 2
    """
}
