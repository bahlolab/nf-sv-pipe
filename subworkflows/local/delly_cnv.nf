
include { DELLY_CNV_CALL     as CALL     } from '../../modules/local/delly_cnv_call'
include { DELLY_MERGE_CNV    as MERGE    } from '../../modules/local/delly_merge_cnv'
include { DELLY_CNV_GENOTYPE as GENOTYPE } from '../../modules/local/delly_cnv_genotype'
include { DELLY_CNV_NORM     as NORM     } from '../../modules/local/delly_cnv_norm'

workflow DELLY_CNV {
    take:
        ref_ch      // value: [ref_fa, [ref_fai, ref_gzi?]]
        fam_bam_ch  // queue: [fam, sam, bam, bai]
        map_ch      // value: [map_fa, map_gzi, map_fai]

    main:
        sam_bam_ch = fam_bam_ch.map { _fam, sam, bam, bai -> [sam, bam, bai] }

        CALL(sam_bam_ch, ref_ch, map_ch)

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

        MERGE(branched.multi.map { fam, _sams, bcfs, csis -> [fam, bcfs, csis] })

        genotype_in = MERGE.out
            .combine(fam_bam_ch.map { fam, sam, bam, bai -> [fam, sam, bam, bai] }, by: 0)
            .map { _fam, sites_bcf, sites_csi, sam, bam, bai -> [sam, bam, bai, sites_bcf, sites_csi] }

        GENOTYPE(genotype_in, ref_ch, map_ch)

        NORM(singleton_bcfs.mix(GENOTYPE.out))

        vcfs = NORM.out.map { sam, bcf, csi -> ['DELLY_CNV', sam, bcf, csi] }

    emit:
        vcfs
}
