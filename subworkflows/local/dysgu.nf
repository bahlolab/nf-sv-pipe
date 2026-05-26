
include { DYSGU_CALL   as CALL   } from '../../modules/local/dysgu_call'
include { DYSGU_TO_BCF as TO_BCF } from '../../modules/local/dysgu_to_bcf'

workflow DYSGU {
    take:
        ref_ch      // value: [ref_fa, [ref_fai, ref_gzi?]]
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        sam_bam_ch = fam_bam_ch.map { _fam, sam, bam, bai -> [sam, bam, bai] }

        CALL(sam_bam_ch, ref_ch)
        TO_BCF(CALL.out)

        vcfs = TO_BCF.out.map { sam, bcf, csi -> ['DYSGU', sam, bcf, csi] }

    emit:
        vcfs
}
