
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

workflow CHORUS {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        chrs_ch     // value channel: List<String> (empty list = no restriction)
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        vcfs = channel.empty()


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

        MATCHA(vcfs, chrs_ch)

    emit:
        vcfs      = vcfs                 // [caller, sam, bcf, csi]
        collapsed = MATCHA.out.collapsed // [sam, bcf, csi]
        merged    = MATCHA.out.merged    // [cohort.bcf, cohort.bcf.csi]
}
