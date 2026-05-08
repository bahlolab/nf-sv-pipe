
process FETCH_REFERENCE_FILES {
    label 'samtools'
    label 'C2M4T4'
    storeDir params.refdir

    output:
    path(smoove_out),                                                                      emit: smoove_excl
    path(delly_excl_out),                                                                  emit: delly_excl
    tuple path(delly_map_out), path("${delly_map_out}.gzi"), path("${delly_map_out}.fai"), emit: delly_map

    script:
    eff            = params.chr_prefix != null ? params.chr_prefix : (params.assembly == 'hg38' ? 'chr' : '')
    smoove_native  = params.assembly == 'hg38' ? 'chr' : ''
    tag            = "${params.assembly}.${eff == 'chr' ? 'chr' : 'nochr'}"
    smoove_out     = "smoove.${tag}.excl.bed"
    delly_excl_out = "delly.${tag}.excl.tsv"
    delly_map_out  = "delly.${tag}.map.gz"

    if (params.assembly == 'hg38') {
        smoove_url     = 'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed'
        delly_excl_url = 'https://raw.githubusercontent.com/dellytools/delly/refs/heads/main/excludeTemplates/human.hg38.excl.tsv'
        delly_map_url  = 'https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz'
    } else {
        smoove_url     = 'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed'
        delly_excl_url = 'https://raw.githubusercontent.com/dellytools/delly/refs/heads/main/excludeTemplates/human.hg19.excl.tsv'
        delly_map_url  = 'https://gear-genomics.embl.de/data/delly/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz'
    }
    """
    # smoove: native has chr iff hg38; strip or add chr only if forced to non-native
    if [ "${eff}" = "${smoove_native}" ]; then
        wget -q -O $smoove_out $smoove_url
    elif [ "${eff}" = "" ]; then
        wget -q -O - $smoove_url | sed 's/^chr//' > $smoove_out
    else
        wget -q -O - $smoove_url | sed '/^#/!s/^/chr/' > $smoove_out
    fi

    # delly excl: always distributed with chr; strip if eff=''
    if [ "${eff}" = "" ]; then
        wget -q -O - $delly_excl_url | sed 's/^chr//' > $delly_excl_out
    else
        wget -q -O $delly_excl_out $delly_excl_url
    fi

    # delly map (bgzip FASTA): always distributed with chr; strip headers if eff=''
    # always generate .gzi and .fai locally
    if [ "${eff}" = "" ]; then
        wget -q -O - $delly_map_url | bgzip -cd  | sed 's/^>chr/>/' | bgzip -@ ${task.cpus} > $delly_map_out
    else
        wget -q -O $delly_map_out $delly_map_url
    fi
    bgzip -r $delly_map_out
    samtools faidx $delly_map_out
    """
}
