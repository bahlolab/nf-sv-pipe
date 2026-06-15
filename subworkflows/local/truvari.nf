
include { TRUVARI_COLLAPSE as COLLAPSE         } from '../../modules/local/truvari_collapse'
include { TRUVARI_MERGE    as MERGE            } from '../../modules/local/truvari_merge'
include { BCF_CLEAN        as CLEAN            } from '../../modules/local/bcf_clean'
include { DUPHOLD                              } from '../../modules/local/duphold'
include { per_sample_by_caller_priority        } from '../../helpers'

workflow TRUVARI {
    take:
        vcfs            // queue: [caller, sam, bcf, csi]
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
                    .map { sam, bcf, csi, bam, bai -> ['TRUVARI', sam, bcf, csi, bam, bai] },
                ref_ch
            )
            duphold = DUPHOLD.out
            fresh   = DUPHOLD.out
        }

        to_merge = fresh.mix(cached_merge_ch)

        MERGE(
            to_merge.map { _sam, bcf, _csi -> bcf }.collect(),
            to_merge.map { _sam, _bcf, csi -> csi }.collect()
        )

        CLEAN(MERGE.out, 'TRUVARI')

    emit:
        duphold                       // [sam, bcf, csi] - freshly duphold-filtered per-sample BCFs (empty when params.duphold=false); feeds the output manifest's TRUVARI-DUPHOLD rows
}
