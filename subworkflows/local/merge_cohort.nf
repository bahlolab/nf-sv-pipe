
include { SVTYPE_SPLIT            as SPLIT         } from '../../modules/local/svtype_split'
include { JASMINE_MERGE           as JASMINE       } from '../../modules/local/jasmine_merge'
include { TRUVARI_COLLAPSE_COHORT as COLLAPSE_COHORT } from '../../modules/local/truvari_collapse_cohort'
include { SV_FINAL_MERGE          as FINAL_MERGE   } from '../../modules/local/sv_final_merge'

workflow MERGE_COHORT {
    take:
        vcfs    // [sam, vcf, tbi]
        ref_ch  // value: [ref_fa, ref_fai]

    main:
        SPLIT(vcfs)

        jasmine_in = SPLIT.out
            .map { sam, jvcf, jtbi, tvcf, ttbi -> [jvcf, jtbi] }
            .collect()
            .map { files ->
                def pairs = files.collate(2)
                [pairs.collect { it[0] }, pairs.collect { it[1] }]
            }

        JASMINE(jasmine_in)

        truvari_in = SPLIT.out
            .map { sam, jvcf, jtbi, tvcf, ttbi -> [tvcf, ttbi] }
            .collect()
            .map { files ->
                def pairs = files.collate(2)
                [pairs.collect { it[0] }, pairs.collect { it[1] }]
            }

        COLLAPSE_COHORT(truvari_in, ref_ch)

        FINAL_MERGE(JASMINE.out.combine(COLLAPSE_COHORT.out))

    emit:
        vcf = FINAL_MERGE.out  // [final_vcf, final_tbi]
}
