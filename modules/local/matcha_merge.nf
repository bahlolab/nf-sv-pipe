
process MATCHA_MERGE {
    label 'bcftools'
    label 'C4M16T4'
    tag   "${chr ?: 'all'}"

    input:
    tuple path(bcfs), path(csis), val(chr)
    val(chr_set)

    output:
    tuple val(chr), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = chr ? "${params.id}.${chr}.cohort.bcf" : "${params.id}.cohort.bcf"
    def chr_arg     = chr     ? "--chrs ${chr}"        : ""
    def chr_set_arg = chr_set ? "--chr-set ${chr_set}" : ""
    """
    export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/local/lib

    matcha merge \\
        --min-jaccard ${params.matcha_min_jaccard} \\
        --threads ${task.cpus} \\
        --missing-to-ref \\
        --write-index \\
        ${chr_arg} \\
        ${chr_set_arg} \\
        -o ${out_bcf} \\
        ${bcfs.join(' ')}
    """
}
