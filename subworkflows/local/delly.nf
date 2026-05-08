
include { DELLY_CALL        as CALL        } from '../../modules/local/delly_call'
include { DELLY_MERGE_SITES as MERGE_SITES } from '../../modules/local/delly_merge_sites'
include { DELLY_GENOTYPE    as GENOTYPE    } from '../../modules/local/delly_genotype'
include { DELLY_BCF_TO_VCF  as BCF_TO_VCF  } from '../../modules/local/delly_bcf_to_vcf'

workflow DELLY {
    take:
        ref_ch      // value: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]
        excl_ch     // value: delly exclude TSV

    main:
        sam_bam_ch = fam_bam_ch.map { _fam, sam, bam, bai -> [sam, bam, bai] }

        CALL(sam_bam_ch, ref_ch, excl_ch)

        call_with_fam = CALL.out
            .join(fam_bam_ch.map { fam, sam, _bam, _bai -> [sam, fam] })
            .map { sam, bcf, csi, fam -> [fam, sam, bcf, csi] }

        fam_calls = call_with_fam.groupTuple(by: 0)

        branched = fam_calls.branch { item ->
            singleton: item[1].size() == 1
            multi:     true
        }

        singleton_bcfs = branched.singleton
            .map { _fam, sams, bcfs, csis -> [sams[0], bcfs[0], csis[0]] }

        MERGE_SITES(branched.multi.map { fam, _sams, bcfs, csis -> [fam, bcfs, csis] })

        genotype_in = MERGE_SITES.out
            .combine(fam_bam_ch.map { fam, sam, bam, bai -> [fam, sam, bam, bai] }, by: 0)
            .map { _fam, sites_bcf, sites_csi, sam, bam, bai -> [sam, bam, bai, sites_bcf, sites_csi] }

        GENOTYPE(genotype_in, ref_ch, excl_ch)

        BCF_TO_VCF(singleton_bcfs.mix(GENOTYPE.out), false)

        vcfs = BCF_TO_VCF.out.map { sam, vcf, csi -> ['DELLY', sam, vcf, csi] }

    emit:
        vcfs
}
