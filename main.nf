#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.ped = ''
params.bams = ''
params.ref_fasta = ''

include { path; read_tsv; get_families; date_ymd } from './nf/functions'
include { manta_call } from './nf/manta_call'
include { annotate_id } from './nf/annotate_id'
include { duphold } from './nf/duphold'
include { bcftools_merge } from './nf/bcftools_merge'
include { filter_duphold } from './nf/filter_duphold'
include { get_pass_ids } from './nf/get_pass_ids'
include { jasmine_merge } from './nf/jasmine_merge'
include { filter_pass_variants } from './nf/filter_pass_variants'

ped = read_tsv(path(params.ped), ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams = read_tsv(path(params.bams), ['iid', 'bam'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')

workflow {

    ref = Channel.value([ref_fa, ref_fai])

    fam_bam_ch =
        Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] } |
        combine(ped.collect { [it.iid, it.fid] }, by: 0) |
        map { it[[3,0,1,2]] }

    fam_vcfs =
        fam_bam_ch |
        map { it[[0,2,3]] } |
        groupTuple(by: 0) |
        combine(ref) |
        manta_call |
        flatMap { it[1].collect { p -> [it[0], p]} } |
        filter { it[1] ==~ '.+diploidSV\\.vcf\\.gz' } |
        annotate_id |
        combine(fam_bam_ch, by:0) |
        map { it[[0,3,1,2,4,5]] } |
        combine(ref) |
        duphold |
        map { it[0,2,3] } |
        groupTuple(by:0) |
        bcftools_merge |
        filter_duphold |
        view

//    fam_vcfs |
//        map { it[1] } |
//        toSortedList() |
//        jasmine_merge |
//        combine(get_pass_ids(fam_vcfs).toSortedList().map { [it] }) |
//        filter_pass_variants
}