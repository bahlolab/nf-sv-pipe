
process BCFTOOLS_MERGE_CALLERS {
    label 'bcftools'
    label 'C2M4T4'
    tag { sam }
    publishDir "${params.progdir}/bcftools_merge_callers", mode: 'symlink'

    input:
        tuple val(sam), val(callers), path(vcfs)

    output:
        tuple val(sam), path(out_vcf), path("${out_vcf}.tbi")

    script:
        out_vcf = "${sam}.callers_merged.vcf.gz"
        """
        # vcfs and callers arrive in priority order (params.callers index)
        callers_arr=($callers)
        vcfs_arr=($vcfs)
        sorted=""
        for i in "\${!vcfs_arr[@]}"; do
            vcf="\${vcfs_arr[\$i]}"
            caller="\${callers_arr[\$i]}"
            name=\$(basename \$vcf)
            echo "\$caller" > sample_name.txt
            bcftools reheader -s sample_name.txt \$vcf -o rh_\${name}
            bcftools sort -Oz --write-index -o sorted_\${name} rh_\${name}
            sorted="\$sorted sorted_\${name}"
            rm rh_\${name}
        done
        bcftools merge --threads $task.cpus -m none \$sorted -Oz -o $out_vcf
        bcftools index --tbi $out_vcf
        rm -f sorted_*.vcf.gz sorted_*.vcf.gz.tbi sample_name.txt
        """
}
