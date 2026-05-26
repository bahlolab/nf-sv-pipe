
process SVDB_COLLAPSE {
    label 'bcftools_svdb'
    label 'C2M4T4'
    tag "$sam"
    publishDir "$params.outdir/svdb_collapse"

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.svdb.collapsed.bcf"
    def priority = callers.join(',')
    def vcf_pairs = [callers, bcfs.collect { it.name.replaceAll(/\.bcf$/, '.vcf.gz') }]
        .transpose()
        .collect { c, v -> "${v}:${c}" }
        .join(' ')
    def filter_cmd = params.svdb_sample_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.svdb_sample_filter}' -Ob -o ${out_bcf} collapsed.bcf"
        : "mv collapsed.bcf ${out_bcf}"
    """
    for BCF in ${bcfs.join(' ')}; do
        out="\${BCF%.bcf}.vcf.gz"
        bcftools view --threads ${task.cpus} -Oz -o "\$out" "\$BCF"
        bcftools index -t --threads ${task.cpus} "\$out"
    done

    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_bnd_distance} \\
        --overlap ${params.svdb_overlap} \\
        --priority ${priority} \\
        --vcf ${vcf_pairs} \\
        | bcftools view --threads ${task.cpus} -Ob -o collapsed.bcf

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
