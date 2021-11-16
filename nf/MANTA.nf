
include { call } from './MANTA/call'
include { annotate_id } from './common/annotate_id'
include { split_sv_types } from './MANTA/split_sv_types'
include { jasmine_merge } from './common/jasmine_merge'
include { get_pass_ids } from './common/get_pass_ids'
include { filter_pass_variants } from './common/filter_pass_variants'
include { bcftools_concat } from './MANTA/bcftools_concat'

//include { duphold } from './MANTA/duphold'
//include { bcftools_merge } from './MANTA/bcftools_merge'
//include { filter_duphold } from './filter_duphold'



workflow MANTA {
    take:
        ref
        fam_bam_ch

    main:

    fam_vcfs = fam_bam_ch |
        map { it[[0,2,3]] } |
        groupTuple(by: 0) |
        take(3) |
        combine(ref) |
        call |
        map { it[0..1] } |
        annotate_id

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
        bcftools_concat

//        combine(fam_bam_ch, by:0) |
//        map { it[[0,3,1,2,4,5]] } |
//        combine(ref) |
//        duphold |
//        map { it[0,2,3] } |
//        groupTuple(by:0) |
//        bcftools_merge |
//        filter_duphold



//    emit:
}
