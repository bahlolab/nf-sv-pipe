
include { manta_call as call } from './manta_call'
include { annotate_id } from './annotate_id'
include { duphold } from './duphold'
include { bcftools_merge } from './bcftools_merge'
include { bcftools_concat } from './bcftools_concat' addParams(pubdir: "output", mode: "copy")
include { filter_duphold } from './filter_duphold'
include { split_sv_types } from './split_sv_types'
include { get_pass_ids } from './get_pass_ids'
include { jasmine_merge } from './jasmine_merge'
include { filter_pass_variants } from './filter_pass_variants'

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
        filter_duphold

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

//    emit:
}
