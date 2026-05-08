
process DELLY_MERGE_SITES {
    label 'delly'
    label 'C2M4T4'
    tag { fam }
    publishDir "${params.progdir}/delly_merge_sites", mode: 'symlink'

    input:
    tuple val(fam), path(bcfs), path(csis)

    output:
    tuple val(fam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${fam}.delly_sites.bcf"
    """
    delly merge -o $out_bcf $bcfs
    """
}
