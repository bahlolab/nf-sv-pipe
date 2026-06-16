
include { SVDB_COLLAPSE   as COLLAPSE         } from '../../modules/local/svdb_collapse'
include { SVDB_MERGE      as MERGE            } from '../../modules/local/svdb_merge'
include { BCFTOOLS_CONCAT as CONCAT           } from '../../modules/local/bcftools_concat'
include { BCF_CLEAN       as CLEAN            } from '../../modules/local/bcf_clean'
include { DUPHOLD                             } from '../../modules/local/duphold'
include { per_sample_by_caller_priority       } from '../../helpers'

workflow SVDB {
    take:
        vcfs            // queue: [caller, sam, bcf, csi]
        chrs_ch         // value channel: List<String> (empty list = no restriction)
        bam_ch          // queue: [sam, bam, bai]
        ref_ch          // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        cached_merge_ch // queue: [sam, bcf, csi] - merge-cache entries to mix into cohort merge (populated in merge-only mode)

    main:
        per_sample = per_sample_by_caller_priority(vcfs)

        COLLAPSE(per_sample)

        duphold = Channel.empty()
        fresh   = COLLAPSE.out
        if (params.duphold) {
            DUPHOLD(
                COLLAPSE.out
                    .combine(bam_ch, by: 0)
                    .map { sam, bcf, csi, bam, bai -> ['SVDB', sam, bcf, csi, bam, bai] },
                ref_ch
            )
            duphold = DUPHOLD.out
            fresh   = DUPHOLD.out
        }

        to_merge = fresh.mix(cached_merge_ch)

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

        sorted_concat_ch = MERGE.out
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { items -> [items.collect { it[1] }, items.collect { it[2] }] }

        concat_split = sorted_concat_ch.multiMap { bcfs, csis ->
            bcfs: bcfs
            csis: csis
        }
        CONCAT(concat_split.bcfs, concat_split.csis, 'SVDB', false)

        CLEAN(CONCAT.out, 'SVDB')

    emit:
        duphold                       // [sam, bcf, csi] - freshly duphold-filtered per-sample BCFs (empty when params.duphold=false); feeds the output manifest's SVDB-DUPHOLD rows
}
