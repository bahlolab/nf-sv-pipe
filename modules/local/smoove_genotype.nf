
process SMOOVE_GENOTYPE {
    label 'smoove'
    label 'C2M8T8'
    tag { sam }

    input:
    tuple val(sam), path(bam), path(bai), path(sites_vcf)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    smoove_vcf = "${sam}-smoove.genotyped.vcf.gz"
    out_bcf    = "${sam}.SMOOVE.bcf"
    """
    smoove genotype $bam -d -x \\
        --processes $task.cpus \\
        --name ${sam} \\
        --outdir . \\
        --fasta $ref_fa \\
        --vcf $sites_vcf
   
    # Soft-filter DEL/DUP by smoove's read-depth fold-change support; other SVTYPEs are untouched.
    bcftools view -i 'SVTYPE="DEL"' $smoove_vcf -Ou | bcftools filter -e 'FMT/DHFFC[0] > 0.7' -s hiDHFFC -Ob -o del.bcf
    bcftools view -i 'SVTYPE="DUP"' $smoove_vcf -Ou | bcftools filter -e 'FMT/DHBFC[0] < 1.3' -s loDHBFC -Ob -o dup.bcf
    bcftools view --threads ${task.cpus} -e 'SVTYPE="DEL" || SVTYPE="DUP"' $smoove_vcf -Ob -o other.bcf 
    ls del.bcf dup.bcf other.bcf | xargs -n1 bcftools index --threads  ${task.cpus}
    bcftools concat -a --threads ${task.cpus} del.bcf dup.bcf other.bcf -Ob -o $out_bcf
    bcftools index --threads ${task.cpus} $out_bcf
    """
}

// If smoove emits duplicate records again, reinstate the dedup workaround
// between `smoove genotype` and `bcftools view`:
//     (
//         bcftools view -h $smoove_vcf
//         bcftools view -H $smoove_vcf | awk '!/^#/ {key = $1$2$3$4$5; if (!seen[key]++) print $0}'
//     ) | bcftools view --threads ${task.cpus} -Ob -o $out_bcf
