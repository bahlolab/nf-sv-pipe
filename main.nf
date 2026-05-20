#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { path               } from './helpers'
include { read_tsv           } from './helpers'
include { check_test_fixtures} from './helpers'
include { check_callers      } from './helpers'
include { get_chrs_ch        } from './helpers'
include { CHORUS             } from './workflows/sv_calling'

workflow {

    check_test_fixtures()
    check_callers()

    def ped     = read_tsv(path(params.ped),  ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
    def bams    = read_tsv(path(params.bams), ['iid', 'bam'])
    def ref_fa  = path(params.ref_fasta)
    def ref_fai = path(params.ref_fasta + '.fai')

    def ref_ch  = Channel.value([ref_fa, ref_fai])
    def chrs_ch = get_chrs_ch()

    def fam_bam_ch =
        Channel.from(bams)
            .map { [it.iid, path(it.bam), path(it.bam + '.bai')] }
            .combine(ped.collect { [it.iid, it.fid] }, by: 0)
            .map { iid, bam, bai, fid -> [params.familial ? fid : iid, iid, bam, bai] }

    CHORUS(ref_ch, chrs_ch, fam_bam_ch)
}
