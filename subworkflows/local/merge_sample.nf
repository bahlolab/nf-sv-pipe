
include { OCTOPUSV_CORRECT        as CORRECT         } from '../../modules/local/octopusv_correct'
include { TRUVARI_COLLAPSE_SAMPLE as COLLAPSE_SAMPLE } from '../../modules/local/truvari_collapse_sample'

workflow MERGE_SAMPLE {
    take:
        vcfs    // queue: [caller, sam, vcf, index]
        ref_ch  // value: [ref_fa, ref_fai]

    main:
        CORRECT(vcfs.map { caller, sam, vcf, _idx -> [caller, sam, vcf] })

        per_sample = CORRECT.out
            .map { caller, sam, vcf -> [sam, caller, vcf] }
            .groupTuple(by: 0)
            .map { sam, callers, vcfs ->
                def order = [callers, vcfs].transpose()
                    .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                def sorted = order.transpose()
                [sam, sorted[0], sorted[1]]
            }

        COLLAPSE_SAMPLE(per_sample, ref_ch)

    emit:
        vcfs = COLLAPSE_SAMPLE.out.map { sam, vcf, _csi -> [sam, vcf] }
}
