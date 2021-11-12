
params.allow_overlap = true
params.pubdir = "progress/bcftools_concat"
params.mode = 'symlink'

process bcftools_concat {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir params.pubdir, mode: params.mode

    input:
        tuple path(vcfs), path(indices)

    output:
        tuple path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = params.id + ".manta-jasmine.vcf.gz"
    """
    bcftools concat $vcfs -Ob -o $out_vcf \\
        ${params.allow_overlap ? '--allow-overlaps' : '--naive-force' }    
    bcftools index -t $out_vcf
    """
}
