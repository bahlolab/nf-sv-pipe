
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
            .map { fam, sam, _bam, _bai -> [fam, sam] }
            .groupTuple(by: 0)
            .map { fam, sams -> [fam, sams.size()] }

        branched = fam_bam_ch
            .combine(fam_sizes, by: 0)
            .branch {
                singleton: it[4] == 1
                multi:     true
            }

        // Singletons: use sam as the call id so output is named per-sample. bam/bai wrapped in lists for MANTA_CALL's bam.join().
        singleton_grouped = branched.singleton
            .map { _fam, sam, bam, bai, _size -> [sam, [bam], [bai]] }

        multi_grouped = branched.multi
            .map { fam, _sam, bam, bai, _size -> [fam, bam, bai] }
            .groupTuple(by: 0)

        id_is_singleton = singleton_grouped.map { sam, _b, _i -> [sam, true] }
            .mix(multi_grouped.map { fam, _b, _i -> [fam, false] })

        CALL(singleton_grouped.mix(multi_grouped), ref_ch, call_regions_ch)

        call_tagged = CALL.out
            .combine(id_is_singleton, by: 0)
            .branch {
                singleton: it[3] == true
                multi:     true
            }

        // Singletons: CALL.out already named per-sample (id=sam) — go straight to FIX.
        singleton_to_fix = call_tagged.singleton
            .map { sam, vcf, tbi, _x -> [sam, vcf, tbi] }

        // Multi: split family-level VCF into per-sample BCFs (intermediate, not published) then FIX.
        SPLIT(call_tagged.multi.map { fam, vcf, tbi, _x -> [fam, vcf, tbi] }, 'MANTA', false)

        split_to_fix = SPLIT.out
            .flatMap { _fam, bcfs, csis ->
                def bs = bcfs instanceof List ? bcfs : [bcfs]
                def cs = csis instanceof List ? csis : [csis]
                [bs, cs].transpose().collect { bcf, csi ->
                    def sam = bcf.name.replaceFirst(/\.MANTA\.bcf$/, '')
                    [sam, bcf, csi]
                }
            }

        FIX_VCF(singleton_to_fix.mix(split_to_fix))

        vcfs = FIX_VCF.out.map { sam, bcf, csi -> ['MANTA', sam, bcf, csi] }

    emit:
        vcfs
}
