
params.caller = 'MANTA'

include { call } from './MANTA/call'
include { convert_inv } from './MANTA/convert_inv'
include { split_sv_types } from './MANTA/split_sv_types'
include { jasmine_merge } from './MANTA/jasmine_merge'
include { duphold } from './MANTA/duphold'
include { vcf_merge } from './MANTA/vcf_merge'
include { filter_duphold } from './common/filter_duphold'
include { publish_vcf } from './common/publish_vcf'
/*
TODO:
    - BND merging, INV and TRA easier (can convert from BND rep)
    - Grapthtyper on merged set?
ISSUES:
    - Unresolved INS not properly merged (ie, with LEFTINSSEQ and RIGHTINSSEQ)
*/

workflow MANTA {
    take:
        ref
        fam_bam_ch

    main:

    manta_vcf = fam_bam_ch |
        map { it[[0,2,3]] } |
        groupTuple(by: 0) |
        combine(ref) |
        call |
//        combine(ref) |
//        convert_inv |
        map { it[1..2] } |
        toSortedList() |
        map { it.transpose() } |
        jasmine_merge |
        combine(fam_bam_ch.map { it.drop(1) }) |
        map { it[2,0,1,3,4] } | //sm, vcf, tbi, bam, bai
        combine(ref) |
        duphold |
        map { it[1].toString() } |
        collectFile(name: 'sample_vcf_list.txt', newLine:true, sort:true) |
        first() |
        vcf_merge |
        map { it[0] } |
        filter_duphold |
        publish_vcf

    emit:
        manta_vcf
}
