
include { TRUVARI_COLLAPSE as COLLAPSE } from '../../modules/local/truvari_collapse'
include { TRUVARI_MERGE    as MERGE    } from '../../modules/local/truvari_merge'
include { DUPHOLD                      } from '../../modules/local/duphold'

workflow TRUVARI {
    take:
        vcfs    // queue: [caller, sam, bcf, csi]
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
        if (params.truvari_duphold) {
            DUPHOLD(COLLAPSE.out.join(bam_ch), ref_ch)
            to_merge = DUPHOLD.out
        }

        MERGE(
            to_merge.map { _sam, bcf, _csi -> bcf }.collect(),
            to_merge.map { _sam, _bcf, csi -> csi }.collect()
        )

    emit:
        collapsed = COLLAPSE.out   // [sam, bcf, csi]
        merged    = MERGE.out      // [cohort.bcf, cohort.bcf.csi]
}
