
include { MATCHA_COLLAPSE as COLLAPSE } from '../../modules/local/matcha_collapse'
include { MATCHA_MERGE    as MERGE    } from '../../modules/local/matcha_merge'

workflow MATCHA {
    take:
        vcfs    // queue: [caller, sam, bcf, csi]

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
        MERGE(
            COLLAPSE.out.map { _sam, bcf, _csi -> bcf }.collect(),
            COLLAPSE.out.map { _sam, _bcf, csi -> csi }.collect()
        )

    emit:
        collapsed = COLLAPSE.out   // [sam, bcf, csi]
        merged    = MERGE.out      // [cohort.bcf, cohort.bcf.csi]
}
