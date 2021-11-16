
include { bins } from './QDNASEQ/bins'
include { mosdepth } from './QDNASEQ/mosdepth'
include { call } from './QDNASEQ/call'
include { annotate_id } from './common/annotate_id' addParams(use_id:false)
include { jasmine_merge } from './common/jasmine_merge' addParams(pubdir: "output", mode: "copy")
include { get_pass_ids } from './common/get_pass_ids'
include { filter_pass_variants } from './common/filter_pass_variants'

workflow QDNASEQ {
    take:
        ref
        fam_bam_ch

    main:
        sam_bam_ch = fam_bam_ch.map { it.drop(1) }

        sample_vcfs = sam_bam_ch |
            take(3) |
            combine(bins()) |
            mosdepth |
            combine(ref.map {it[1]} ) |
            call |
            annotate_id

        sample_vcfs |
            map { it[1] } |
            toSortedList() |
            map { ['QDNASEQ', it] } |
            jasmine_merge |
            map { it[1] } |
            combine(get_pass_ids(sample_vcfs).toSortedList().map { [it] }) |
            first() |
            filter_pass_variants

//    emit:
}
