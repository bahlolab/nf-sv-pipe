
process vep {
    cpus 4
    memory '4 GB'
    time '1 h'
    publishDir "output/vep", mode: 'copy'

    input:
        tuple path(vcf), path(fasta), path(fai), path(cache)

    output:
        path(vep_vcf)

    script:
    vep_vcf = vcf.name.replaceAll('.vcf.gz', '.vep.vcf.gz')
    vep_output_opts = [
        '--ccds',
        '--symbol',
        '--numbers',
        '--protein',
        '--variant_class'
    ].join(' ')
    vep_filter_opts = [
        '--pick_allele_gene'
    ].join(' ')
    """
    bcftools view --no-version  $vcf |
        vep --input_file STDIN \\
            $vep_output_opts \\
            $vep_filter_opts \\
            --fork 4 \\
            --format vcf \\
            --vcf \\
            --cache \\
            --offline \\
            --no_stats \\
            --fasta $fasta \\
            --assembly $params.vep_assembly \\
            --cache_version $params.vep_cache_ver \\
            --dir $cache \\
            --output_file STDOUT |
            bcftools view --no-version -Oz -o $vep_vcf 
    """
}
