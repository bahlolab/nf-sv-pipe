
process jasmine {
    cpus 4
    memory '4 GB'
    time '1 h'
    publishDir "progress/jasmine", mode: 'symlink'

    input:
    path(vcfs)

    output:
    path(out_vcf)

    script:
    out_vcf = params.id + '.merged.vcf.gz'
    """
    gzip -fdk *diploidSV.id.vcf.gz
    echo *diploidSV.id.vcf | tr ' ' '\\n' > vcf_list.txt
    
    jasmine \\
        file_list=vcf_list.txt \\
        out_file=merged.vcf \\
        threads=$task.cpus \\
        --output_genotypes 
        # --normalize_type # causes vcf format errors with bcftools
    
    bcftools sort merged.vcf -Oz -o sorted.vcf.gz
    bcftools query -l sorted.vcf.gz | sed 's:^\\([0-9]*_\\)\\(.*\\):\\1\\2 \\2:' > rename.txt
    bcftools reheader sorted.vcf.gz --sample rename.txt --output $out_vcf
    
    rm merged.vcf sorted.vcf.gz *diploidSV.id.vcf output -r
    """
}

