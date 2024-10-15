
process get_exclude_bed {
    cpus 1
    memory '1 GB'
    time '1 h'
    publishDir "${params.progdir}/SMOOVE/exclude_bed", mode: 'symlink'
    container null

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
