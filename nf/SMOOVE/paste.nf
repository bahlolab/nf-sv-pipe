
process paste {
    cpus 1
    memory '4 GB'
    time '1 h'
    publishDir 'progress/SMOOVE/paste', mode: 'symlink'
    container 'quay.io/biocontainers/smoove:0.2.5--0'

    input:
        tuple path(vcfs), path(indices)

    output:
        path(vcf)

    script:
        vcf = "${params.id}.smoove.square.vcf.gz"
        """
        mkdir tmp && export TMPDIR=tmp
        smoove paste $vcfs --name $params.id
        mv $vcf tmp.vcf.gz
        bcftools +fill-tags tmp.vcf.gz -Oz -o $vcf -- -t AF,AC,AN
        """
}
