#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id        = ''
params.ped       = ''
params.bams      = ''
params.ref_fasta = ''
params.assembly  = 'hg38'
// Order of params.callers also defines caller priority for per-sample merging:
// the first caller's records are preferred by truvari --keep maxqual.
params.callers   = ['MANTA', 'SMOOVE', 'CNVNATOR']
params.copy_ref  = false
params.copy_bams = false
params.outdir    = 'output'
params.progdir   = 'progress'
params.bin_size  = 1000
params.chrs      = null

params.truvari_intra_refdist  = 500
params.truvari_intra_pctseq   = 0.7
params.truvari_intra_pctsize  = 0.7
params.truvari_intra_bnddist  = 500
params.jasmine_max_dist       = 500
params.truvari_cohort_refdist = 500
params.truvari_cohort_pctseq  = 0.7
params.truvari_cohort_pctsize = 0.7
params.truvari_cohort_bnddist = 500

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
include { SV_CALLING     } from './workflows/sv_calling'

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
            .map { iid, bam, bai, fid -> [fid, iid, bam, bai] }

    SV_CALLING(ref_ch, fam_bam_ch)
}
