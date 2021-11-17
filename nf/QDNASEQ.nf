
params.caller = 'QDNASEQ'

include { bins } from './QDNASEQ/bins'
include { mosdepth } from './QDNASEQ/mosdepth'
include { call } from './QDNASEQ/call'
include { set_id } from './common/set_id' addParams(use_id:false)
include { jasmine_merge } from './common/jasmine_merge' addParams(pubdir: "output", mode: "copy")
include { get_pass_ids } from './common/get_pass_ids'
include { filter_pass_variants } from './common/filter_pass_variants'
include { publish_vcf } from './common/publish_vcf'

workflow QDNASEQ {
    take:
        ref
        fam_bam_ch

    main:
        sam_bam_ch = fam_bam_ch.map { it.drop(1) }

        sample_vcfs = sam_bam_ch |
            combine(bins()) |
            mosdepth |
            combine(ref.map {it[1]} ) |
            call |
            set_id

        qdnaseq_vcf = sample_vcfs |
            map { it[1] } |
            toSortedList() |
            map { ['ALL', it] } |
            jasmine_merge |
            map { it[1] } |
            combine(get_pass_ids(sample_vcfs).toSortedList().map { [it] }) |
            first() |
            filter_pass_variants |
            publish_vcf

    emit:
        qdnaseq_vcf
}
