
process JASMINE_MERGE {
    label 'jasminesv'
    label 'C4M32T4'
    publishDir "${params.progdir}/jasmine_merge", mode: 'symlink'

    input:
    path(vcfs)

    output:
    path(out_vcf)

    script:
    out_vcf = "${params.id}.cohort_jasmine.vcf.gz"
    """
    tmp_vcfs=""
    for vcf in $vcfs; do
        name=\$(basename \$vcf .vcf.gz)
        gzip -dc \$vcf > \${name}.vcf
        echo "\${name}.vcf" >> filelist.txt
        tmp_vcfs="\$tmp_vcfs \${name}.vcf"
    done
    jasmine file_list=filelist.txt out_file=cohort_jasmine.vcf \\
        threads=$task.cpus \\
        max_dist=$params.jasmine_max_dist \\
        --output_genotypes \\
        --keep_var_ids
    rm -f \$tmp_vcfs filelist.txt
    gzip -c cohort_jasmine.vcf > $out_vcf
    rm cohort_jasmine.vcf
    """
}
