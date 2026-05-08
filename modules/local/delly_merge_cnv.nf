
process DELLY_MERGE_CNV {
    label 'delly'
    label 'C2M4T4'
    tag { fam }
    publishDir "${params.progdir}/delly_merge_cnv", mode: 'symlink'

    input:
    tuple val(fam), path(bcfs), path(csis)

    output:
    tuple val(fam), path(out_bcf), path("${out_bcf}.csi")

    script:
    out_bcf = "${fam}.delly_cnv_sites.bcf"
    """
    delly merge -e -p -m 1000 -n 100000 -o $out_bcf $bcfs
    """
}
