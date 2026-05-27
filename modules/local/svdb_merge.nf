
process SVDB_MERGE {
    label 'bcftools_svdb'
    label 'C2M4T24'
    tag "${chr ?: 'all'}"

    input:
    tuple path(bcfs), path(csis), val(chr)

    output:
    tuple val(chr), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = chr ? "${params.id}.${chr}.SVDB.merge.bcf" : "${params.id}.SVDB.merge.bcf"
    def region_arg = chr ? "-r ${chr}" : ""
    def filter_cmd = params.svdb_cohort_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.svdb_cohort_filter}' -Ob -o ${out_bcf} merged.bcf"
        : "mv merged.bcf ${out_bcf}"
    """
    for BCF in ${bcfs.join(' ')}; do
        out="\${BCF%.bcf}.vcf.gz"
        bcftools view --threads ${task.cpus} ${region_arg} -Oz -o "\$out" "\$BCF"
        bcftools index -t --threads ${task.cpus} "\$out"
    done

    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_bnd_distance} \\
        --overlap ${params.svdb_overlap} \\
        --vcf ${bcfs.collect { it.name.replaceAll(/\.bcf$/, '.vcf.gz') }.join(' ')} \\
        | bcftools view --threads ${task.cpus} -Ob -o merged.bcf

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
