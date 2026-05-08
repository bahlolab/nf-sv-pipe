#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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

include { path; read_tsv } from './nf/functions'
include { CHORUS         } from './workflows/sv_calling'

ped    = read_tsv(path(params.ped),       ['fid', 'iid', 'pid', 'mid', 'sex', 'phe'])
bams   = read_tsv(path(params.bams),      ['iid', 'bam'])
ref_fa  = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')

workflow {

    ref_ch = Channel.value([ref_fa, ref_fai])

    fam_bam_ch =
        Channel.from(bams)
            .map { [it.iid, path(it.bam), path(it.bam + '.bai')] }
            .combine(ped.collect { [it.iid, it.fid] }, by: 0)
            .map { iid, bam, bai, fid -> [params.familial ? fid : iid, iid, bam, bai] }

    CHORUS(ref_ch, fam_bam_ch)
}
