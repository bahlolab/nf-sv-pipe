
include { get_exclude_bed } from './SMOOVE/exclude_bed'
include { call } from './SMOOVE/call'
include { merge } from './SMOOVE/merge'
include { genotype } from './SMOOVE/genotype'
include { paste } from './SMOOVE/paste'
include { filter_fc } from './SMOOVE/filter_fc'


workflow SMOOVE {

    take:
        ref
        fam_bam_ch

    main:
        sam_bam_ch = fam_bam_ch.map { it.drop(1) }

        smoove_vcf = sam_bam_ch |
            combine(ref) |
            combine(get_exclude_bed()) |
//            filter { ['S33366_1', 'S34233_2', 'S33365_1'].contains(it[0]) } |
            call |
            map { it.drop(1) } |
            toSortedList() |
            map { it.transpose() } |
            combine(ref) |
            first() |
            merge |
            combine(ref) |
            combine(sam_bam_ch) |
//            combine(sam_bam_ch.filter{['S33366_1', 'S34233_2', 'S33365_1'].contains(it[0])}) |
            genotype |
            map { it.drop(1) } |
            toSortedList() |
            map { it.transpose() } |
            paste |
            filter_fc

    emit:
        smoove_vcf
}
