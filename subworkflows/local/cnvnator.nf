
include { CNVNATOR_PROCESS_REF as PROCESS_REF } from '../../modules/local/cnvnator_process_ref'
include { CNVNATOR_CALL        as CALL        } from '../../modules/local/cnvnator_call'
include { CNVNATOR_TO_BCF      as TO_BCF      } from '../../modules/local/cnvnator_to_bcf'

workflow CNVNATOR {
    take:
        ref_ch      // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        chrs_ch     // value channel: List<String> (empty list = no restriction)
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        sam_bam_ch = fam_bam_ch.map { fam, sam, bam, bai -> [sam, bam, bai] }

        PROCESS_REF(
            ref_ch,
            chrs_ch
        )

        proc_ref = PROCESS_REF.out

        chrs_str_ch = chrs_ch.map { it ? it.join(' ') : '' }

        CALL(
            sam_bam_ch,
            proc_ref.map { ref_dir, fai -> ref_dir },
            chrs_str_ch
        )

        TO_BCF(CALL.out, proc_ref.map { ref_dir, fai -> fai })

        vcfs = TO_BCF.out
            .map { sam, bcf, csi -> ['CNVNATOR', sam, bcf, csi] }

    emit:
        vcfs
}
