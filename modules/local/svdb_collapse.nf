
process SVDB_COLLAPSE {
    label 'bcftools_svdb'
    label 'C2M4T4'
    tag "$sam"
    publishDir "$params.outdir/SVDB"

    input:
    tuple val(sam), val(callers), path(bcfs), path(csis)

    output:
    tuple val(sam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${sam}.SVDB.bcf"
    def priority = callers.join(',')
    def bases = [callers, bcfs.collect { it.name.replaceAll(/\.bcf$/, '') }].transpose()
    // NB: subset suffixes use '_itvl'/'_bnd', NOT '.itvl'/'.bnd'. svdb ignores the VCF spec and
    // embeds characters from the input filename into the INFO field IDs it adds; a '.' in the name
    // produces invalid INFO IDs and breaks downstream bcftools. Keep these to [A-Za-z0-9_] only.
    def itvl_pairs = bases.collect { c, b -> "${b}_itvl.vcf.gz:${c}" }.join(' ')
    def bnd_pairs  = bases.collect { c, b -> "${b}_bnd.vcf.gz:${c}" }.join(' ')
    def filter_cmd = params.svdb_sample_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.svdb_sample_filter}' -Ob -o ${out_bcf} collapsed.bcf && rm -f collapsed.bcf"
        : "mv collapsed.bcf ${out_bcf}"
    """
    # split each per-caller VCF into interval (DEL/DUP/INV) and BND/INS subsets.
    # suffix with '_itvl'/'_bnd' (underscore, not '.'): svdb embeds filename chars into the
    # INFO IDs it adds, and a '.' yields invalid INFO IDs that break downstream bcftools.
    for BCF in ${bcfs.join(' ')}; do
        base="\${BCF%.bcf}"
        bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \\
            --threads ${task.cpus} -Oz -o "\${base}_itvl.vcf.gz" "\$BCF"
        bcftools index -t --threads ${task.cpus} "\${base}_itvl.vcf.gz"
        bcftools view -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="INS"' \\
            --threads ${task.cpus} -Oz -o "\${base}_bnd.vcf.gz" "\$BCF"
        bcftools index -t --threads ${task.cpus} "\${base}_bnd.vcf.gz"
    done

    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_itvl_bnd_distance} \\
        --overlap ${params.svdb_itvl_overlap} \\
        --priority ${priority} \\
        --vcf ${itvl_pairs} \\
        | bcftools view --threads ${task.cpus} -Ob -o itvl.bcf
    bcftools index --threads ${task.cpus} itvl.bcf

    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_bnd_bnd_distance} \\
        --overlap ${params.svdb_bnd_overlap} \\
        --priority ${priority} \\
        --vcf ${bnd_pairs} \\
        | bcftools view --threads ${task.cpus} -Ob -o bnd.bcf
    bcftools index --threads ${task.cpus} bnd.bcf

    rm -f *_itvl.vcf.gz *_itvl.vcf.gz.tbi *_bnd.vcf.gz *_bnd.vcf.gz.tbi

    bcftools concat --allow-overlaps --threads ${task.cpus} \\
        itvl.bcf bnd.bcf -Ob -o collapsed.bcf
    rm -f itvl.bcf itvl.bcf.csi bnd.bcf bnd.bcf.csi

    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
