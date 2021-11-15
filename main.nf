#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.ped = ''
params.bams = ''
params.ref_fasta = ''
params.assembly = 'hg38'

include { path; read_tsv; get_families; date_ymd } from './nf/functions'
include { QDNASEQ } from './nf/QDNASEQ'
include { MANTA } from './nf/MANTA'

ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')

workflow {

    ref_ch = Channel.value([ref_fa, ref_fai])

    fam_bam_ch =
        Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] }

//    MANTA(ref_ch, fam_bam_ch)
    QDNASEQ(ref_ch, fam_bam_ch)

}