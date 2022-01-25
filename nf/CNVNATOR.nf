
params.caller = 'CNVNATOR'
params.bin_size = 1000
params.chrs = null

include { process_ref } from './CNVNATOR/process_ref'
include { call } from './CNVNATOR/call'
include { to_vcf } from './CNVNATOR/to_vcf'
include { jasmine_merge } from './common/jasmine_merge'
include { publish_vcf } from './common/publish_vcf'
/*
    TODO - sample filtering for samples that have high noise and too many calls
 */
workflow CNVNATOR {
    take:
        ref
        fam_bam_ch

    main:
        chrs = params.chrs ?:
            params.assembly == 'hg38' ?
                ((1..22) + ['X', 'Y']).collect{ 'chr' + it } :
                ((1..22) + ['X', 'Y']).collect{ is.toString() }

        sam_bam_ch = fam_bam_ch.map { it.drop(1) }

        mod_ref = ref |
            map { it + [chrs] } |
            process_ref

        cnvnator_vcf = sam_bam_ch |
            map { it + [chrs.join(' ')] } |
            combine(process_ref.out.map{ it[0] }) |
            call |
            combine(process_ref.out) |
            to_vcf |
            map { it[1..2] } |
            toSortedList() |
            map { it.transpose() } |
            jasmine_merge |
            publish_vcf

    emit:
        cnvnator_vcf
}
