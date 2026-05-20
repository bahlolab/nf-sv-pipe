
include { SMOOVE_CALL     as CALL     } from '../../modules/local/smoove_call'
include { SMOOVE_MERGE    as MERGE    } from '../../modules/local/smoove_merge'
include { SMOOVE_GENOTYPE as GENOTYPE } from '../../modules/local/smoove_genotype'

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

        // Singletons: use the sample's own CALL output as the sites VCF for GENOTYPE
        singleton_sites = branched.singleton
            .map { fam, sams, vcfs, csis -> [fam, vcfs[0]] }

        // Multi-member families: merge per-sample calls into a family sites VCF
        MERGE(
            branched.multi.map { fam, sams, vcfs, csis -> [fam, vcfs, csis] },
            ref_ch
        )

        // Unified per-family sites VCF channel
        sites_ch = singleton_sites.mix(MERGE.out)

        // Pair each sample with its family's sites VCF and re-genotype
        genotype_in = sites_ch
            .combine(fam_bam_ch.map { fam, sam, bam, bai -> [fam, sam, bam, bai] }, by: 0)
            .map { fam, sites_vcf, sam, bam, bai -> [sam, bam, bai, sites_vcf] }

        GENOTYPE(genotype_in, ref_ch)

        vcfs = GENOTYPE.out.map { sam, bcf, csi -> ['SMOOVE', sam, bcf, csi] }

    emit:
        vcfs
}
