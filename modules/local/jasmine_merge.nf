
process JASMINE_MERGE {
    label 'jasminesv'
    label 'C4M32T4'
    publishDir "${params.progdir}/jasmine_merge", mode: 'symlink'

    input:
        tuple path(vcfs), path(tbis)

    output:
        tuple path(out_vcf), path("${out_vcf}.csi")

    script:
        out_vcf = "${params.id}.cohort_jasmine.vcf.gz"
        """
        ls $vcfs > filelist.txt
        jasmine file_list=filelist.txt out_file=cohort_jasmine.vcf \\
            threads=$task.cpus \\
            max_dist=$params.jasmine_max_dist \\
            --output_genotypes \\
            --keep_var_ids
        bgzip cohort_jasmine.vcf
        bcftools index cohort_jasmine.vcf.gz
        mv cohort_jasmine.vcf.gz $out_vcf
        mv cohort_jasmine.vcf.gz.csi ${out_vcf}.csi
        """
}
