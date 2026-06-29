
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

    // svdb --merge requires VCF.gz + TBI; it also can't handle sample names with dots/slashes,
    // so BCFs are renamed to SAM_<safe> on conversion. TBI is created on the _itvl/_bnd subsets
    // (which svdb actually reads), not on the intermediate SAM_*.vcf.gz (only bcftools view reads that).
    // svdb adds set/FOUNDBY/svdb_origin/SUPP_VEC/VARID during per-sample collapse (SVDB_COLLAPSE);
    // stripping them here prevents duplicate INFO tags when svdb appends them again at cohort level.
    def svdb_tags = 'INFO/set,INFO/FOUNDBY,INFO/svdb_origin,INFO/SUPP_VEC,INFO/VARID'
    def safe_bases = bcfs.collect { bcf ->
        "SAM_" + bcf.name.replaceAll(/\..*$/, '').replaceAll(/[^A-Za-z0-9_]/, '_')
    }
    // NB: subset suffixes use '_itvl'/'_bnd', NOT '.itvl'/'.bnd'. svdb ignores the VCF spec and
    // embeds characters from the input filename into the INFO field IDs it adds; a '.' in the name
    // produces invalid INFO IDs (e.g. 'svdb_origin.itvl') and breaks downstream bcftools. Keep these
    // to [A-Za-z0-9_] only.
    def itvl_vcfs = safe_bases.collect { it + '_itvl.vcf.gz' }.join(' ')
    def bnd_vcfs  = safe_bases.collect { it + '_bnd.vcf.gz' }.join(' ')

    // strip INFO bloat inline in each svdb merge pipe, before bcftools concat — invalid/duplicate
    // INFO fields (N*M per-sample entries) would otherwise break concat's header union step
    def info_strip_pipe = params.svdb_info_keep
        ? "| bcftools annotate -x '^${params.svdb_info_keep.collect { 'INFO/' + it }.join(',')}' -Ob -o"
        : "| bcftools view --threads ${task.cpus} -Ob -o"

    // optional PASS/quality filter on the cohort-merged result
    def filter_cmd = params.svdb_cohort_filter
        ? "bcftools view --threads ${task.cpus} -i '${params.svdb_cohort_filter}' -Ob -o ${out_bcf} merged.bcf && rm -f merged.bcf"
        : "mv merged.bcf ${out_bcf}"
    """
    # convert, strip stale SVDB tags, and split into SVTYPE subsets per sample.
    # safe name: strip extension (%%.*), replace non-[A-Za-z0-9_] with '_', prefix SAM_ —
    # must match the Groovy safe_bases computation so itvl_vcfs/bnd_vcfs resolve correctly.
    # suffix with '_itvl'/'_bnd' (underscore, not '.'): svdb embeds filename chars into the
    # INFO IDs it adds, and a '.' yields invalid INFO IDs that break downstream bcftools.
    for BCF in ${bcfs.join(' ')}; do
        stem="\${BCF%%.*}"
        safe="SAM_\${stem//[^A-Za-z0-9_]/_}"
        bcftools annotate ${region_arg} -x ${svdb_tags} -Ou "\$BCF" \\
            | bcftools view --threads ${task.cpus} -Oz -o "\${safe}.vcf.gz"
        bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV"' \\
            --threads ${task.cpus} -Oz -o "\${safe}_itvl.vcf.gz" "\${safe}.vcf.gz"
        bcftools index -t --threads ${task.cpus} "\${safe}_itvl.vcf.gz"
        bcftools view -i 'INFO/SVTYPE="BND" || INFO/SVTYPE="INS"' \\
            --threads ${task.cpus} -Oz -o "\${safe}_bnd.vcf.gz" "\${safe}.vcf.gz"
        bcftools index -t --threads ${task.cpus} "\${safe}_bnd.vcf.gz"
        rm -f "\${safe}.vcf.gz"
    done

    # cohort merge — interval subset; strip INFO bloat inline before writing BCF
    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_itvl_bnd_distance} \\
        --overlap ${params.svdb_itvl_overlap} \\
        --vcf ${itvl_vcfs} \\
        ${info_strip_pipe} itvl.bcf
    bcftools index --threads ${task.cpus} itvl.bcf

    # cohort merge — BND/INS subset; strip INFO bloat inline before writing BCF
    svdb \\
        --merge \\
        --bnd_distance ${params.svdb_bnd_bnd_distance} \\
        --overlap ${params.svdb_bnd_overlap} \\
        --vcf ${bnd_vcfs} \\
        ${info_strip_pipe} bnd.bcf
    bcftools index --threads ${task.cpus} bnd.bcf

    rm -f SAM_*_itvl.vcf.gz SAM_*_itvl.vcf.gz.tbi SAM_*_bnd.vcf.gz SAM_*_bnd.vcf.gz.tbi

    bcftools concat --allow-overlaps --threads ${task.cpus} \\
        itvl.bcf bnd.bcf -Ob -o merged.bcf

    # optional cohort filter
    ${filter_cmd}
    bcftools index --threads ${task.cpus} ${out_bcf}
    """
}
