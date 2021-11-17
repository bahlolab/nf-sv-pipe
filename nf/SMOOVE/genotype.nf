
process genotype {
    cpus 2
    memory '4 GB'
    time '4 h'
    tag { sam }
    publishDir 'progress/SMOOVE/genotype', mode: 'symlink'
    container 'quay.io/biocontainers/smoove:0.2.5--0'

    input:
        tuple path(vcf), path(ref_fa), path(ref_fai), val(sam), path(bam), path(bai)

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
            --vcf $vcf 
        """
}
