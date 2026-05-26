
include { MANTA                 } from '../subworkflows/local/manta'
include { SMOOVE                } from '../subworkflows/local/smoove'
include { CNVNATOR              } from '../subworkflows/local/cnvnator'
include { DELLY                 } from '../subworkflows/local/delly'
include { DELLY_CNV             } from '../subworkflows/local/delly_cnv'
include { DYSGU                 } from '../subworkflows/local/dysgu'
include { FETCH_REFERENCE_FILES } from '../modules/local/fetch_reference_files'
include { MAKE_CALL_REGIONS     } from '../modules/local/make_call_regions'
include { COPY_BAMS             } from '../modules/local/copy_bams'
include { PASS_FILTER           } from '../modules/local/pass_filter'
include { MATCHA                } from '../subworkflows/local/matcha'
include { TRUVARI               } from '../subworkflows/local/truvari'

workflow SVPLEX {
    take:
        ref_ch      // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        chrs_ch     // value channel: List<String> (empty list = no restriction)
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        vcfs = channel.empty()

        fam_bam_ch = fam_bam_ch
            .groupTuple(by:0)
            .map {fam, sam, bam, bai -> [ groupKey(fam, sam.size()), sam, bam, bai] }
            .transpose()

        def needs_ref = params.callers.intersect(['MANTA', 'SMOOVE', 'DELLY', 'DELLY_CNV'])
        if (needs_ref) { FETCH_REFERENCE_FILES() }

        if (params.copy_bams) { fam_bam_ch = COPY_BAMS(fam_bam_ch) }

        if (params.callers.contains('MANTA')) {
            call_regions_ch = MAKE_CALL_REGIONS(ref_ch, chrs_ch, FETCH_REFERENCE_FILES.out.delly_excl)
            vcfs = vcfs.mix(MANTA(ref_ch, fam_bam_ch, call_regions_ch).vcfs)
        }
        if (params.callers.contains('SMOOVE'))    { vcfs = vcfs.mix(SMOOVE(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.smoove_excl).vcfs) }
        if (params.callers.contains('CNVNATOR'))  { vcfs = vcfs.mix(CNVNATOR(ref_ch, chrs_ch, fam_bam_ch).vcfs) }
        if (params.callers.contains('DELLY'))     { vcfs = vcfs.mix(DELLY(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.delly_excl).vcfs) }
        if (params.callers.contains('DELLY_CNV')) { vcfs = vcfs.mix(DELLY_CNV(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.delly_map).vcfs) }
        if (params.callers.contains('DYSGU'))     { vcfs = vcfs.mix(DYSGU(ref_ch, fam_bam_ch).vcfs) }

        branched = vcfs.branch { caller, sam, bcf, csi ->
            filter:    params.apply_filters.contains(caller)
            no_filter: true
        }
        vcfs = PASS_FILTER(branched.filter).mix(branched.no_filter)

        matcha_collapsed  = channel.empty()
        matcha_merged     = channel.empty()
        truvari_collapsed = channel.empty()
        truvari_merged    = channel.empty()

        sam_bam_ch = fam_bam_ch.map { _fam, sam, bam, bai -> [sam, bam, bai] }

        if (params.matcha) {
            MATCHA(vcfs, chrs_ch, sam_bam_ch, ref_ch)
            matcha_collapsed = MATCHA.out.collapsed
            matcha_merged    = MATCHA.out.merged
        }
        if (params.truvari) {
            TRUVARI(vcfs, sam_bam_ch, ref_ch)
            truvari_collapsed = TRUVARI.out.collapsed
            truvari_merged    = TRUVARI.out.merged
        }

    emit:
        vcfs              = vcfs              // [caller, sam, bcf, csi]
        matcha_collapsed  = matcha_collapsed  // [sam, bcf, csi]
        matcha_merged     = matcha_merged     // [cohort.bcf, cohort.bcf.csi]
        truvari_collapsed = truvari_collapsed // [sam, bcf, csi]
        truvari_merged    = truvari_merged    // [cohort.bcf, cohort.bcf.csi]
}
