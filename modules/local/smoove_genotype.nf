
process SMOOVE_GENOTYPE {
    label 'smoove'
    label 'C2M8T2'
    tag { sam }
    publishDir "${params.progdir}/SMOOVE/genotype", mode: 'symlink'

    input:
    tuple val(sam), path(bam), path(bai), path(sites_vcf)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path(out_vcf), path("${out_vcf}.csi")

    script:
    out_vcf = "${sam}-smoove.genotyped.vcf.gz"
    """
    mkdir tmp && export TMPDIR=tmp
    smoove genotype $bam -d -x \\
        --processes $task.cpus \\
        --name ${sam} \\
        --outdir . \\
        --fasta $ref_fa \\
        --vcf $sites_vcf
    # work around for smoove duplicate record bug
    (
        bcftools view -h $out_vcf
        bcftools view -H $out_vcf | awk '!/^#/ {key = \$1\$2\$3\$4\$5; if (!seen[key]++) print \$0}'
    ) | bcftools view -Oz -o tmp.vcf.gz
    mv tmp.vcf.gz $out_vcf
    bcftools index -f $out_vcf
    """
}
