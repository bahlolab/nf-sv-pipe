
include { SVTYPE_SPLIT            as SPLIT               } from '../../modules/local/svtype_split'
include { JASMINE_MERGE           as JASMINE             } from '../../modules/local/jasmine_merge'
include { BCFTOOLS_SORT_INDEX     as SORT_INDEX_JASMINE  } from '../../modules/local/bcftools_sort_index'
include { BCFTOOLS_MERGE_COHORT   as MERGE_TRUVARI_SAMPLES } from '../../modules/local/bcftools_merge_cohort'
include { TRUVARI_COLLAPSE_COHORT as COLLAPSE_COHORT     } from '../../modules/local/truvari_collapse_cohort'
include { SV_FINAL_MERGE          as FINAL_MERGE         } from '../../modules/local/sv_final_merge'

workflow MERGE_COHORT {
    take:
        vcfs    // [sam, vcf]
        ref_ch  // value: [ref_fa, ref_fai]

    main:
        SPLIT(vcfs)

        jasmine_in = SPLIT.out
            .map { sam, jvcf, tvcf -> jvcf }
            .collect()

        JASMINE(jasmine_in)
        SORT_INDEX_JASMINE(JASMINE.out)

        truvari_in = SPLIT.out
            .map { sam, jvcf, tvcf -> tvcf }
            .collect()

        MERGE_TRUVARI_SAMPLES(truvari_in)
        COLLAPSE_COHORT(MERGE_TRUVARI_SAMPLES.out, ref_ch)

        FINAL_MERGE(SORT_INDEX_JASMINE.out.combine(COLLAPSE_COHORT.out))

    emit:
        vcf = FINAL_MERGE.out  // [final_vcf, final_tbi]
}
