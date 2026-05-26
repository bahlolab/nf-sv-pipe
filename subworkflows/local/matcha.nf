
include { MATCHA_COLLAPSE as COLLAPSE } from '../../modules/local/matcha_collapse'
include { MATCHA_MERGE    as MERGE    } from '../../modules/local/matcha_merge'
include { DUPHOLD                     } from '../../modules/local/duphold'

workflow MATCHA {
    take:
        vcfs    // queue: [caller, sam, bcf, csi]
        chrs_ch // value channel: List<String> (empty list = no restriction)
        bam_ch  // queue: [sam, bam, bai]
        ref_ch  // value channel: [ref_fa, [ref_fai, ref_gzi?]]

    main:
        chrs_str_ch = chrs_ch.map { it ? it.join(',') : '' }

        per_sample = vcfs
            .map { caller, sam, bcf, csi -> [groupKey(sam.toString(), params.callers.size()), caller, bcf, csi] }
            .groupTuple(by:0)
            .map { sam, callers, bcfs, csis ->
                def order = [callers, bcfs, csis].transpose()
                    .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                def sorted = order.transpose()
                [sam.target, sorted[0], sorted[1], sorted[2]]
            }

        COLLAPSE(per_sample, chrs_str_ch)

        to_merge = COLLAPSE.out
        if (params.matcha_duphold) {
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
