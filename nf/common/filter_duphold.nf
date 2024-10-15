
params.min_len = 500
params.max_del_fc = 0.70
params.min_dup_fc = 1.25

process filter_duphold {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "${params.progdir}/${params.caller}/filter_duphold", mode: 'symlink'

    input:
        path(vcf)

    output:
        tuple path(out_vcf), path("${out_vcf}.tbi")

    script:

    out_vcf = vcf.name.replaceAll('.vcf.gz', '.filter-duphold.vcf.gz')
    """
    filter_duphold_fc.py \\
        --input $vcf \\
        --output soft-filterd.vcf.gz \\
        --max-del-fc $params.max_del_fc \\
        --min-dup-fc $params.min_dup_fc \\
        --min-len $params.min_len
    bcftools view -f PASS soft-filterd.vcf.gz -Oz -o $out_vcf
    bcftools index -t $out_vcf
    """
}