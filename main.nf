#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.ped = ''
params.bams = ''
params.ref_fasta = ''
params.assembly = 'hg38'
params.callers = ['MANTA']
//params.callers = ['MANTA', 'QDNASEQ', 'SMOOVE', 'CNVNATOR']
params.copy_ref = true
params.copy_bams = true

include { path; read_tsv; get_families; date_ymd } from './nf/functions'
include { copy_ref } from './nf/common/copy_ref'
include { copy_bams } from './nf/common/copy_bams'
include { MANTA } from './nf/MANTA'
include { QDNASEQ } from './nf/QDNASEQ'
include { SMOOVE } from './nf/SMOOVE'
include { CNVNATOR } from './nf/CNVNATOR'

ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')

workflow {

    ref_ch = Channel.value([ref_fa, ref_fai])
    if (params.copy_ref) { ref_ch = copy_ref(ref_ch) }

    fam_bam_ch =
        Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] }
    if (params.copy_bams) { fam_bam_ch = copy_bams(fam_bam_ch) }

    if (params.callers.contains('SMOOVE')) {
        SMOOVE(ref_ch, fam_bam_ch)
    }
    if (params.callers.contains('MANTA')) {
        MANTA(ref_ch, fam_bam_ch)
    }
    if (params.callers.contains('QDNASEQ')) {
        QDNASEQ(ref_ch, fam_bam_ch)
    }
    if (params.callers.contains('CNVNATOR')) {
        CNVNATOR(ref_ch, fam_bam_ch)
    }
}