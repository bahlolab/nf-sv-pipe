
include { MANTA        } from '../subworkflows/local/manta'
include { SMOOVE       } from '../subworkflows/local/smoove'
include { CNVNATOR     } from '../subworkflows/local/cnvnator'
include { MERGE_SAMPLE } from '../subworkflows/local/merge_sample'
include { MERGE_COHORT } from '../subworkflows/local/merge_cohort'

workflow SV_CALLING {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        vcfs = Channel.empty()

        if (params.callers.contains('MANTA'))    { vcfs = vcfs.mix(MANTA(ref_ch, fam_bam_ch).vcfs) }
        if (params.callers.contains('SMOOVE'))   { vcfs = vcfs.mix(SMOOVE(ref_ch, fam_bam_ch).vcfs) }
        if (params.callers.contains('CNVNATOR')) { vcfs = vcfs.mix(CNVNATOR(ref_ch, fam_bam_ch).vcfs) }

        MERGE_SAMPLE(vcfs, ref_ch)
        MERGE_COHORT(MERGE_SAMPLE.out.vcfs, ref_ch)

    emit:
        vcfs           = vcfs                    // [caller, sample_id, vcf, vcf_index]
        merged_sample  = MERGE_SAMPLE.out.vcfs   // [sam, vcf]
        merged_cohort  = MERGE_COHORT.out.vcf    // [final_vcf, final_tbi]
}
