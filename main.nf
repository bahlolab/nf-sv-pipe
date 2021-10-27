#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
Overview:
    1. call SVs per family with manta
    2. merge SVs across cohort with Jasmine
    3. Annotate SVs with VEP?
 */
//imput params
params.id = ''
params.ped = ''
params.bams = ''
// ref params
params.ref_fasta = ''
//params.vep_cache = ''
//params.vep_cache_ver = ''
//params.vep_assembly = params.ref_hg38 ? 'GRCh38' : 'GRCh37'
// exec params

include { path; read_tsv; get_families; date_ymd } from './nf/functions'

ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')
//vep_cache = path(params.vep_cache)

workflow {

    ref = Channel.value([ref_fa, ref_fai])

    Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] } |
        groupTuple(by: 0) |
        combine(ref) |
        view

}