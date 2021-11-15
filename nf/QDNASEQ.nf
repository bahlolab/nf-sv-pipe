
include { qdnaseq_bins } from './qdnaseq_bins'
include { mosdepth } from './mosdepth'
include { qdnaseq_call as call } from './qdnaseq_call'
include { annotate_id } from './annotate_id' addParams(use_id:false)
include { jasmine_merge } from './jasmine_merge' addParams(pubdir: "output", mode: "copy")


workflow QDNASEQ {
    take:
        ref
        fam_bam_ch

    main:
        fam_bam_ch |
            combine(qdnaseq_bins()) |
            mosdepth |
            combine(ref.map {it[1]} ) |
            call |
            map { it[1..2] } |
            annotate_id |
            map { it[1] } |
            toSortedList() |
            map { ['QDNASEQ', it] } |
            jasmine_merge |
            view

//    emit:
}
