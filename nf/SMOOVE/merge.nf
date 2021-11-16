
process merge {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir 'progress/smoove_merge', mode: 'symlink'
    container null
    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/smoove'

    input:
        tuple path(vcfs), path(indices), path(ref_fa), path(ref_fai)

    output:
        path(vcf)

    script:
        pref = "${params.id}.smooth-merged"
        vcf = "${pref}.sites.vcf.gz"
        """
        mkdir tmp && export TMPDIR=tmp
        smoove merge $vcfs \\
            --name $pref \\
            --fasta $ref_fa
        """
}
