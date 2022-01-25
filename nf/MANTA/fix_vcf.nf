/*
    Fix SVLEN should be length 1 in header
 */
process fix_vcf {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/$params.caller/fix_vcf", mode: 'symlink'
    tag { id }

    input:
    tuple val(id), path(vcf), path(tbi)

    output:
    tuple val(id), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = vcf.name.replaceAll('.vcf.gz', '.fixed.vcf.gz')
    """
    bcftools view -h $vcf | 
        sed 's:##INFO=<ID=SVLEN,Number=.,:##INFO=<ID=SVLEN,Number=1,:' > header.txt
    bcftools reheader $vcf -h header.txt -o $out_vcf
    bcftools index -t $out_vcf
    """
}

