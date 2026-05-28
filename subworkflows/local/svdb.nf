
include { SVDB_COLLAPSE   as COLLAPSE } from '../../modules/local/svdb_collapse'
include { SVDB_MERGE      as MERGE    } from '../../modules/local/svdb_merge'
include { BCFTOOLS_CONCAT as CONCAT   } from '../../modules/local/bcftools_concat'
include { BCF_CLEAN       as CLEAN    } from '../../modules/local/bcf_clean'
include { DUPHOLD                     } from '../../modules/local/duphold'

workflow SVDB {
    take:
        vcfs    // queue: [caller, sam, bcf, csi]
        chrs_ch // value channel: List<String> (empty list = no restriction)
        bam_ch  // queue: [sam, bam, bai]
        ref_ch  // value channel: [ref_fa, [ref_fai, ref_gzi?]]

    main:
        per_sample = vcfs
            .map { caller, sam, bcf, csi -> [groupKey(sam.toString(), params.callers.size()), caller, bcf, csi] }
            .groupTuple(by:0)
            .map { sam, callers, bcfs, csis ->
                def order = [callers, bcfs, csis].transpose()
                    .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                def sorted = order.transpose()
                [sam.target, sorted[0], sorted[1], sorted[2]]
            }

        COLLAPSE(per_sample)

        to_merge = COLLAPSE.out
        if (params.duphold) {
            DUPHOLD(
                COLLAPSE.out
                    .combine(bam_ch, by: 0)
                    .map { sam, bcf, csi, bam, bai -> ['SVDB', sam, bcf, csi, bam, bai] },
                ref_ch
            )
            to_merge = DUPHOLD.out
        }

        merge_input = to_merge
            .map { sm, bcf, csi -> [true, sm, bcf, csi] }
            .groupTuple(by: 0)
            .map { _key, sms, bcfs, csis ->
                def sorted = [sms, bcfs, csis].transpose().sort { a, b -> a[0] <=> b[0] }
                def s = sorted.transpose()
                [s[1], s[2]]
            }
            .combine(chrs_ch.map { it ?: [null] }.flatten())

        MERGE(merge_input)

        sorted_concat_ch = MERGE.out   // [chr, bcf, csi]
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { items -> [items.collect { it[1] }, items.collect { it[2] }] }

        concat_split = sorted_concat_ch.multiMap { bcfs, csis ->
            bcfs: bcfs
            csis: csis
        }
        CONCAT(concat_split.bcfs, concat_split.csis, 'SVDB')

        CLEAN(CONCAT.out, 'SVDB', params.svdb_info_keep)

    emit:
        collapsed = COLLAPSE.out   // [sam, bcf, csi]
        merged    = CLEAN.out      // [cohort.bcf, cohort.bcf.csi]
}
