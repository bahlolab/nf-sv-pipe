
include { OCTOPUSV_CORRECT       as CORRECT       } from '../../modules/local/octopusv_correct'
include { BCFTOOLS_MERGE_CALLERS as MERGE_CALLERS } from '../../modules/local/bcftools_merge_callers'
include { TRUVARI_COLLAPSE_SAMPLE as COLLAPSE_SAMPLE } from '../../modules/local/truvari_collapse_sample'

workflow MERGE_SAMPLE {
    take:
        vcfs    // [caller, sam, vcf, index]
        ref_ch  // value: [ref_fa, ref_fai]

    main:
        CORRECT(vcfs.map { caller, sam, vcf, idx -> [caller, sam, vcf] })

        per_sample = CORRECT.out
            .map { caller, sam, vcf -> [sam, caller, vcf] }
            .groupTuple(by: 0)
            .map { sam, callers, vcfs ->
                def order = [callers, vcfs].transpose()
                    .sort { a, b -> params.callers.indexOf(a[0]) <=> params.callers.indexOf(b[0]) }
                def sorted = order.transpose()
                [sam, sorted[0], sorted[1]]
            }

        MERGE_CALLERS(per_sample)

        COLLAPSE_SAMPLE(MERGE_CALLERS.out, ref_ch)

    emit:
        vcfs = COLLAPSE_SAMPLE.out  // [sam, vcf]
}
