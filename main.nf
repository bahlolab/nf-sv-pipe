#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { path               } from './helpers'
include { read_tsv           } from './helpers'
include { check_test_fixtures} from './helpers'
include { check_callers      } from './helpers'
include { check_apply_filters} from './helpers'
include { get_chrs_ch        } from './helpers'
include { SVPLEX             } from './workflows/sv_plex'

workflow {

    check_test_fixtures()
    check_callers()
    check_apply_filters()

    def ped     = read_tsv(path(params.ped),  ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
    def bams    = read_tsv(path(params.bams), ['iid', 'bam'])
    def ref_fa  = path(params.ref_fasta)
    def ref_idx = [path(params.ref_fasta + '.fai')]
    if (params.ref_fasta.endsWith('.gz')) {
        ref_idx << path(params.ref_fasta + '.gzi')
    }

    def ref_ch  = Channel.value([ref_fa, ref_idx])
    def chrs_ch = get_chrs_ch()

    def fam_bam_ch =
        Channel.from(bams)
            .map { [it.iid, path(it.bam), path(it.bam + '.bai')] }
            .combine(ped.collect { [it.iid, it.fid] }, by: 0)
            .map { iid, bam, bai, fid -> [params.familial ? fid : iid, iid, bam, bai] }

    SVPLEX(ref_ch, chrs_ch, fam_bam_ch)
}
