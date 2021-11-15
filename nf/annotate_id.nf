params.use_id = true

process annotate_id {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "progress/annotate_id", mode: 'symlink'
    tag {id}

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path(out_vcf), path("${out_vcf}.tbi")

    script:
    out_vcf = vcf.name.replaceAll('.vcf.gz', '.id.vcf.gz')
    """
    bcftools annotate $vcf -Oz -o $out_vcf \\
        --set-id '${id}_${params.use_id ? '%ID' : '%SVTYPE\\_%POS\\_%END' }'
    bcftools index -t $out_vcf
    """
}

