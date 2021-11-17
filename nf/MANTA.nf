
params.caller = 'MANTA'

include { call } from './MANTA/call'
include { set_id } from './common/set_id'
include { split_sv_types } from './MANTA/split_sv_types'
include { jasmine_merge } from './common/jasmine_merge'
include { get_pass_ids } from './common/get_pass_ids'
include { filter_pass_variants } from './common/filter_pass_variants'
include { vcf_concat } from './MANTA/vcf_concat'
include { duphold } from './MANTA/duphold'
include { vcf_merge } from './MANTA/vcf_merge'
include { filter_duphold } from './common/filter_duphold'
include { publish_vcf } from './common/publish_vcf'

workflow MANTA {
    take:
        ref
        fam_bam_ch

    main:

    fam_vcfs = fam_bam_ch |
        map { it[[0,2,3]] } |
        groupTuple(by: 0) |
        combine(ref) |
        call |
        map { it[0..1] } |
        set_id

    fam_vcfs |
        split_sv_types |
        flatMap { it[1] instanceof List ? it[1] : [it[1]] } |
        map { [(it.name =~ /([A-Z]+)\.vcf\.gz$/)[0][1], it] } |
        groupTuple(by:0) |
        jasmine_merge |
        map { it[1] } |
        combine(get_pass_ids(fam_vcfs).toSortedList().map { [it] }) |
        filter_pass_variants |
        toList() |
        map { it.transpose() } |
        vcf_concat |
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

//    emit:
}
