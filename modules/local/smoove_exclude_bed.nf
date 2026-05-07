
process SMOOVE_EXCLUDE_BED {
    label 'smoove'
    label 'C2M2T2'
    publishDir "${params.progdir}/SMOOVE/exclude_bed", mode: 'symlink'

    output:
    path(bed)

    script:
    url = params.assembly == 'hg38' ?
        'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed' :
        'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed'
    bed = url.replaceAll('.+/', '')
    """
    wget $url
    """
}
