
process SMOOVE_CALL {
    label 'smoove'
    label 'C2M16T4'
    tag { sam }
    publishDir "${params.progdir}/SMOOVE/call", mode: 'symlink'

    input:
        tuple val(sam), path(bam), path(bai)
        tuple path(ref_fa), path(ref_fai)
        path(exclude)

    output:
        tuple val(sam), path(vcf), path("${vcf}.csi")

    script:
        vcf = "${sam}-smoove.genotyped.vcf.gz"
        """
        mkdir tmp && export TMPDIR=tmp
        /usr/bin/time -v smoove call $bam \\
            --processes $task.cpus \\
            --outdir . \\
            --exclude $exclude \\
            --name $sam \\
            --fasta $ref_fa \\
            --genotype
        """
}
