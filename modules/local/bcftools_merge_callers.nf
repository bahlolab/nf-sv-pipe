
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
    # caller[0] keeps the original sample name; subsequent callers are renamed to caller name
    set -- $callers
    sorted=""
    for vcf in $vcfs; do
        caller=\$1; shift
        name=\$(basename \$vcf)
        bcftools sort -Oz -o sorted_\${name} \$vcf
        echo "\$caller" > sample_name.txt
        bcftools reheader -s sample_name.txt sorted_\${name} -o rh_\${name}
        bcftools index --tbi rh_\${name}
        sorted="\$sorted rh_\${name}"
        rm sorted_\${name}
    done
    bcftools merge --threads $task.cpus -m none \$sorted -Oz -o $out_vcf
    bcftools index --tbi $out_vcf
    rm -f sorted_*.vcf.gz sorted_*.vcf.gz.tbi rh_*.vcf.gz rh_*.vcf.gz.tbi sample_name.txt
    """
}
