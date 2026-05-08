
process TRUVARI_COLLAPSE_SAMPLE {
    label 'truvari'
    label 'C2M4T4'
    tag { sam }
    publishDir "${params.progdir}/truvari_collapse_sample", mode: 'symlink'

    input:
    tuple val(sam), val(callers), path(vcfs)
    tuple path(ref_fa), path(ref_fai)

    output:
    tuple val(sam), path("${sam}.consensus.vcf.gz"), path("${sam}.consensus.vcf.gz.tbi")

    script:
    def hires_set   = ['MANTA', 'DELLY', 'SMOOVE'] as Set
    def pairs       = [callers, vcfs].transpose()
    def hires_args  = pairs.findAll { name, __ -> hires_set.contains(name) }
        .collect { name, vcf -> "${name}:${vcf}" }.join(' ')
    def lowres_args = pairs.findAll { name, __ -> !hires_set.contains(name) }
        .collect { name, vcf -> "${name}:${vcf}" }.join(' ')
    def lowres_flag = lowres_args ? "--low-res ${lowres_args}" : ""
    """
    truvari_collapse_staged.py \\
        --ref         ${ref_fa} \\
        --output      ${sam}.consensus.vcf.gz \\
        --high-res    ${hires_args} \\
        ${lowres_flag} \\
        --hires-args  "${params.truvari_hires_args}" \\
        --lowres-args "${params.truvari_lowres_args}"
    """
}
