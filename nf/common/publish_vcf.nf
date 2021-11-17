
process publish_vcf {
    cpus 1
    memory '1 GB'
    time '1 m'
    container null
    executor 'local'
    publishDir "output", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    tuple path(out_vcf), path(out_tbi)

    script:
    out_vcf = "${params.id}.${params.caller}.vcf.gz"
    out_tbi = out_vcf + '.tbi'
    """
    ln -s `readlink $vcf` $out_vcf
    ln -s `readlink $tbi` $out_tbi
    """
}

