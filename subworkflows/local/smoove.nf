
include { SMOOVE_CALL        as CALL        } from '../../modules/local/smoove_call'
include { SMOOVE_MERGE       as MERGE       } from '../../modules/local/smoove_merge'
include { SMOOVE_GENOTYPE    as GENOTYPE    } from '../../modules/local/smoove_genotype'
include { SMOOVE_TO_BCF      as TO_BCF      } from '../../modules/local/smoove_to_bcf'

workflow SMOOVE {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]
        excl_ch     // value: smoove exclude BED

    main:
        sam_bam_ch = fam_bam_ch.map { fam, sam, bam, bai -> [sam, bam, bai] }

        CALL(sam_bam_ch, ref_ch, excl_ch)

        // Annotate call output with family ID via join on sample key
        call_with_fam = CALL.out
            .join(fam_bam_ch.map { fam, sam, bam, bai -> [sam, fam] })
            .map { sam, vcf, csi, fam -> [fam, sam, vcf, csi] }

        // Group per family, then branch by family size
        fam_calls = call_with_fam.groupTuple(by: 0)

        branched = fam_calls.branch {
            singleton: it[1].size() == 1
            multi:     true
        }

        // Singletons: use CALL output directly (already genotyped via --genotype)
        singleton_vcfs = branched.singleton
            .map { fam, sams, vcfs, csis -> [sams[0], vcfs[0], csis[0]] }

        // Multi-member families: one merge per family, then re-genotype each sample
        MERGE(
            branched.multi.map { fam, sams, vcfs, csis -> [fam, vcfs, csis] },
            ref_ch
        )

        // Pair each sample in multi-member families with their family's sites VCF
        genotype_in = MERGE.out
            .combine(fam_bam_ch.map { fam, sam, bam, bai -> [fam, sam, bam, bai] }, by: 0)
            .map { fam, sites_vcf, sam, bam, bai -> [sam, bam, bai, sites_vcf] }

        GENOTYPE(genotype_in, ref_ch)

        TO_BCF(singleton_vcfs.mix(GENOTYPE.out))

        vcfs = TO_BCF.out.map { sam, bcf, csi -> ['SMOOVE', sam, bcf, csi] }

    emit:
        vcfs
}
