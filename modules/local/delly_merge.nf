
process DELLY_MERGE {
    label 'delly'
    label 'C2M4T4'
    tag { fam }

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
