
include { MANTA_CALL         as CALL         } from '../../modules/local/manta_call'
include { MANTA_FIX_VCF      as FIX_VCF      } from '../../modules/local/manta_fix_vcf'
include { MANTA_SPLIT_SAMPLE as SPLIT_SAMPLE } from '../../modules/local/manta_split_sample'

workflow MANTA {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        fam_sizes = fam_bam_ch
            .map { fam, sam, bam, bai -> [fam, sam] }
            .groupTuple(by: 0)
            .map { fam, sams -> [fam, sams.size()] }

        branched = fam_bam_ch
            .combine(fam_sizes, by: 0)
            .branch {
                singleton: it[4] == 1
                multi:     true
            }

        // Singletons: use sam as the call id so output is named after the sample.
        // Wrap bam/bai in lists so bam.join() in MANTA_CALL works uniformly.
        singleton_grouped = branched.singleton
            .map { fam, sam, bam, bai, size -> [sam, [bam], [bai]] }

        multi_grouped = branched.multi
            .map { fam, sam, bam, bai, size -> [fam, bam, bai] }
            .groupTuple(by: 0)

        // Tag each call id as singleton (true) or multi (false) for downstream routing
        id_is_singleton = singleton_grouped.map { sam, bams, bais -> [sam, true] }
            .mix(multi_grouped.map    { fam, bams, bais -> [fam, false] })

        CALL(singleton_grouped.mix(multi_grouped), ref_ch)
        FIX_VCF(CALL.out)

        fix_tagged = FIX_VCF.out
            .join(id_is_singleton)
            .branch {
                singleton: it[3] == true
                multi:     true
            }

        // Singletons: id IS sam and file is already named per-sample — emit directly
        singleton_vcfs = fix_tagged.singleton
            .map { sam, vcf, tbi, _x -> ['MANTA', sam, vcf, tbi] }

        // Multi: pair each family VCF with every sample in that family, then split
        split_in = fix_tagged.multi
            .map { fam, vcf, tbi, _x -> [fam, vcf, tbi] }
            .combine(fam_bam_ch.map { fam, sam, bam, bai -> [fam, sam] }, by: 0)

        SPLIT_SAMPLE(split_in)

        vcfs = singleton_vcfs.mix(SPLIT_SAMPLE.out.map { sam, vcf, tbi -> ['MANTA', sam, vcf, tbi] })

    emit:
        vcfs
}
