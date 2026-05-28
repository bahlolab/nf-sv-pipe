
process PUBLISH {
    publishDir { "${params.caller_vcf_dir}/${caller}" }, mode: 'copy'

    input:
    tuple val(caller), val(sam), path(bcf), path(csi)

    script:
    'true'
}
