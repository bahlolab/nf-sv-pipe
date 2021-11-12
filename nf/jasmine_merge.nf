
process jasmine_merge {
    cpus 2
    memory '4 GB'
    time '1 h'
    publishDir "progress/jasmine_merge", mode: 'symlink'
    tag {type}

    input:
    tuple val(type), path(vcfs)

    output:
    tuple val(type), path(out_vcf)

    script:
    out_vcf = params.id + '.' + type + '.merged.vcf.gz'
    """
    gzip -fdk $vcfs
    echo $vcfs | sed 's:\\.gz::g' | tr ' ' '\\n' > vcf_list.txt
    
    jasmine \\
        file_list=vcf_list.txt \\
        out_file=merged.vcf \\
        threads=$task.cpus \\
        --output_genotypes \\
        --default_zero_genotype
        # --normalize_type # causes vcf format errors with bcftools
    
    bcftools sort merged.vcf -Ou |
        bcftools +fill-tags -Oz -o sorted.vcf.gz
    bcftools query -l sorted.vcf.gz | sed 's:^\\([0-9]*_\\)\\(.*\\):\\1\\2 \\2:' > rename.txt
    bcftools reheader sorted.vcf.gz --sample rename.txt --output $out_vcf
    
    rm sorted.vcf.gz *.vcf output -r
    """
}

