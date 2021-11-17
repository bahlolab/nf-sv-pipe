
params.caller = 'SMOOVE'

include { get_exclude_bed } from './SMOOVE/exclude_bed'
include { call } from './SMOOVE/call'
include { merge } from './SMOOVE/merge'
include { genotype } from './SMOOVE/genotype'
include { paste } from './SMOOVE/paste'
include { filter_duphold } from './common/filter_duphold'
include { publish_vcf } from './common/publish_vcf'

workflow SMOOVE {

    take:
        ref
        fam_bam_ch

    main:
        sam_bam_ch = fam_bam_ch.map { it.drop(1) }

        smoove_vcf = sam_bam_ch |
            combine(ref) |
            combine(get_exclude_bed()) |
            call |
            map { it.drop(1) } |
            toSortedList() |
            map { it.transpose() } |
            combine(ref) |
            first() |
            merge |
            combine(ref) |
            combine(sam_bam_ch) |
            genotype |
            map { it.drop(1) } |
            toSortedList() |
            map { it.transpose() } |
            paste |
            filter_duphold |
            publish_vcf

    emit:
        smoove_vcf
}
