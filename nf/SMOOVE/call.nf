
process call {
    cpus 2
    memory '8 GB'
    time '6 h'
    tag { sam }
    publishDir 'progress/smoove_call', mode: 'symlink'
    container null
    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/smoove'
//    container 'quay.io/biocontainers/smoove:0.2.5--0'

    input:
        tuple val(sam), path(bam), path(bai), path(ref_fa), path(ref_fai), path(exclude)

    output:
        tuple val(sam), path(vcf), path("${vcf}.csi")

    script:
        vcf = "${sam}-smoove.genotyped.vcf.gz"
        """
        mkdir tmp && export TMPDIR=tmp
        smoove call $bam \\
            --processes $task.cpus \\
            --outdir . \\
            --exclude $exclude \\
            --name $sam \\
            --fasta $ref_fa \\
            --genotype
        """
}
