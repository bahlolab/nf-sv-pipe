
include { DELLY_CALL     as CALL     } from '../../modules/local/delly_call'
include { DELLY_MERGE    as MERGE    } from '../../modules/local/delly_merge'
include { DELLY_GENOTYPE as GENOTYPE } from '../../modules/local/delly_genotype'
include { DELLY_NORM     as NORM     } from '../../modules/local/delly_norm'

workflow DELLY {
    take:
        ref_ch      // value: [ref_fa, [ref_fai, ref_gzi?]]
        fam_bam_ch  // queue: [fam, sam, bam, bai]
        excl_ch     // value: delly exclude TSV

    main:
        fam_sizes = fam_bam_ch
            .map { fam, sam, _bam, _bai -> [fam, sam] }
            .groupTuple(by: 0)
            .map { fam, sams -> [fam, sams.size()] }

        sam_bam_ch = fam_bam_ch
            .combine(fam_sizes, by: 0)
            .map { _fam, sam, bam, bai, _size -> [sam, bam, bai] }

        CALL(sam_bam_ch, ref_ch, excl_ch)

        call_with_fam = CALL.out
            .combine(fam_bam_ch.map { fam, sam, _bam, _bai -> [sam, fam] }, by: 0)
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

        GENOTYPE(genotype_in, ref_ch, excl_ch)

        NORM(singleton_bcfs.mix(GENOTYPE.out))

        vcfs = NORM.out.map { sam, bcf, csi -> ['DELLY', sam, bcf, csi] }

    emit:
        vcfs
}
