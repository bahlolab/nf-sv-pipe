
process filter_pass_variants {
    label 'pysam'
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/filter_pass_variants", mode: 'symlink'

    input:
        tuple path(vcf), path(id_files)

    output:
        path(filtered)

    script:
    filtered = vcf.name.replaceAll('.vcf.gz', '.filtered.vcf.gz')
    """
    bcftools view -Ou $vcf |
        filter_pass_variants.py \\
            --ids ${id_files.join(' --ids ') } \\
            --out $filtered
    """
}

