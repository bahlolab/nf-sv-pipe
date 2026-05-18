
include { MATCHA_COLLAPSE } from '../../modules/local/matcha_collapse'
include { MATCHA_MERGE    } from '../../modules/local/matcha_merge'

workflow MERGE {
    take:
        vcfs    // queue: [caller, sam, bcf, csi]

    main:
        per_sample = vcfs
            .map { caller, sam, bcf, csi -> [sam, caller, bcf, csi] }
            .groupTuple(by: 0)
            .map { sam, callers, bcfs, csis ->
                def order = [callers, bcfs, csis].transpose()
                    .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                def sorted = order.transpose()
                [sam, sorted[0], sorted[1], sorted[2]]
            }

        MATCHA_COLLAPSE(per_sample)
        MATCHA_MERGE(
            MATCHA_COLLAPSE.out.map { sam, bcf, csi -> bcf }.collect(),
            MATCHA_COLLAPSE.out.map { sam, bcf, csi -> csi }.collect()
        )

    emit:
        collapsed = MATCHA_COLLAPSE.out   // [sam, bcf, csi]
        merged    = MATCHA_MERGE.out      // [cohort.bcf, cohort.bcf.csi]
}
