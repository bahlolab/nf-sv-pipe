
process paste {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir 'progress/smoove_paste', mode: 'symlink'
    container null
    conda '/stornext/Home/data/allstaff/m/munro.j/miniconda3/envs/smoove'

    input:
        tuple path(vcfs), path(indices)

    output:
        path(vcf)

    script:
        vcf = "${params.id}.smoove.square.vcf.gz"
        """
        mkdir tmp && export TMPDIR=tmp
        smoove paste $vcfs --name $params.id
        """
}
