params.pubdir = "progress/jasmine_merge"
params.mode = 'symlink'

process jasmine_merge {
    cpus 2
    memory '4 GB'
    time '1 h'
    publishDir params.pubdir, mode: params.mode
    tag { uid }

    input:
    tuple val(uid), path(vcfs)

    output:
    tuple val(uid), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = params.id + '.' + uid + '.merged.vcf.gz'
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
    bcftools index -t $out_vcf
    
    rm sorted.vcf.gz *.vcf output -r
    """
}

