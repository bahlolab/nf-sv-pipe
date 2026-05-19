#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { path; read_tsv } from './nf/functions'
include { CHORUS         } from './workflows/sv_calling'

workflow {

    if (workflow.profile.contains('test')) {
        [params.bams, params.ped, params.ref_fasta].each { f ->
            if (!file(f).exists()) {
                error """\
                Test fixture missing: ${f}
                Generate fixtures first by running:
                    bash test/generate_fixtures.sh
                """.stripIndent()
            }
        }
    }

    def supported_callers = ['MANTA', 'SMOOVE', 'CNVNATOR', 'DELLY', 'DELLY_CNV']
    def unsupported = params.callers - supported_callers
    if (unsupported) {
        error "Unsupported callers in params.callers: ${unsupported}. Supported: ${supported_callers}"
    }
    def dups = params.callers.countBy { it }.findAll { k, n -> n > 1 }.keySet()
    if (dups) {
        error "Duplicate callers in params.callers: ${dups}"
    }

    def ped    = read_tsv(path(params.ped),       ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
    def bams   = read_tsv(path(params.bams),      ['iid', 'bam'])
    def ref_fa  = path(params.ref_fasta)
    def ref_fai = path(params.ref_fasta + '.fai')

    def ref_ch = Channel.value([ref_fa, ref_fai])

    def fam_bam_ch =
        Channel.from(bams)
            .map { [it.iid, path(it.bam), path(it.bam + '.bai')] }
            .combine(ped.collect { [it.iid, it.fid] }, by: 0)
            .map { iid, bam, bai, fid -> [params.familial ? fid : iid, iid, bam, bai] }

    CHORUS(ref_ch, fam_bam_ch)
}
