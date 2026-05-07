
include { CNVNATOR_PROCESS_REF as PROCESS_REF } from '../../modules/local/cnvnator_process_ref'
include { CNVNATOR_CALL        as CALL        } from '../../modules/local/cnvnator_call'
include { CNVNATOR_TO_VCF      as TO_VCF      } from '../../modules/local/cnvnator_to_vcf'

workflow CNVNATOR {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        chrs = params.chrs ?:
            (params.assembly == 'hg38'
                ? ((1..22) + ['X', 'Y']).collect { 'chr' + it }
                : ((1..22) + ['X', 'Y']).collect { it.toString() })

        sam_bam_ch = fam_bam_ch.map { fam, sam, bam, bai -> [sam, bam, bai] }

        PROCESS_REF(ref_ch.map { fa, fai -> [fa, fai, chrs] })

        proc_ref = PROCESS_REF.out

        CALL(
            sam_bam_ch.map { sam, bam, bai -> [sam, bam, bai, chrs.join(' ')] },
            proc_ref.map { ref_dir, fai -> ref_dir }
        )

        TO_VCF(CALL.out, proc_ref)

        vcfs = TO_VCF.out
            .map { sam, vcf, csi -> ['CNVNATOR', sam, vcf, csi] }

    emit:
        vcfs
}
