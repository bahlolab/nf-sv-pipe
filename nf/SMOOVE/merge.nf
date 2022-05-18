
process merge {
    cpus 1
    memory '4 GB'
    time '1 h'
    publishDir 'progress/SMOOVE/merge', mode: 'symlink'
    container 'quay.io/biocontainers/smoove:0.2.5--0'

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
