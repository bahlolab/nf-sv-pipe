
include { SMOOVE_CALL    as CALL  } from '../../modules/local/smoove_call'
include { BCFTOOLS_SPLIT as SPLIT } from '../../modules/local/bcftools_split'

workflow SMOOVE {
    take:
        ref_ch      // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        fam_bam_ch  // queue: [fam, sam, bam, bai]
        excl_ch     // value: smoove exclude BED

    main:
        fam_bams_ch = fam_bam_ch
            .map { fam, sam, bam, bai -> [fam, bam, bai] }
            .groupTuple()

        CALL(fam_bams_ch, ref_ch, excl_ch)

        SPLIT(CALL.out, 'SMOOVE', true)

        vcfs = SPLIT.out
            .flatMap { fam, bcfs, csis ->
                def bs = bcfs instanceof List ? bcfs : [bcfs]
                def cs = csis instanceof List ? csis : [csis]
                [bs, cs].transpose().collect { bcf, csi ->
                    def sam = bcf.name.replaceFirst(/\.SMOOVE\.bcf$/, '')
                    ['SMOOVE', sam, bcf, csi]
                }
            }

    emit:
        vcfs
}
