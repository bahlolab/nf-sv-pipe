
include { MANTA_CALL     as CALL    } from '../../modules/local/manta_call'
include { MANTA_FIX_VCF  as FIX_VCF } from '../../modules/local/manta_fix_vcf'
include { BCFTOOLS_SPLIT as SPLIT   } from '../../modules/local/bcftools_split'

workflow MANTA {
    take:
        ref_ch          // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        fam_bam_ch      // queue: [fam, sam, bam, bai]
        call_regions_ch // value channel: [call_regions.bed.gz, call_regions.bed.gz.tbi]

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

        CALL(singleton_grouped.mix(multi_grouped), ref_ch, call_regions_ch)

        FIX_VCF(CALL.out)

        fix_tagged = FIX_VCF.out
            .combine(id_is_singleton, by: 0)
            .branch {
                singleton: it[3] == true
                multi:     true
            }

        // Singletons: id IS sam and file is already named per-sample — emit directly
        singleton_vcfs = fix_tagged.singleton
            .map { sam, bcf, csi, _x -> ['MANTA', sam, bcf, csi] }

        // Multi: split family-level BCF into per-sample BCFs
        SPLIT(fix_tagged.multi.map { fam, bcf, csi, _x -> [fam, bcf, csi] }, 'MANTA')

        split_vcfs = SPLIT.out
            .flatMap { fam, bcfs, csis ->
                def bs = bcfs instanceof List ? bcfs : [bcfs]
                def cs = csis instanceof List ? csis : [csis]
                [bs, cs].transpose().collect { bcf, csi ->
                    def sam = bcf.name.replaceFirst(/\.MANTA\.bcf$/, '')
                    ['MANTA', sam, bcf, csi]
                }
            }

        vcfs = singleton_vcfs.mix(split_vcfs)

    emit:
        vcfs
}
