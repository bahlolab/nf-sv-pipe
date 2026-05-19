
include { MANTA                  } from '../subworkflows/local/manta'
include { SMOOVE                 } from '../subworkflows/local/smoove'
include { CNVNATOR               } from '../subworkflows/local/cnvnator'
include { DELLY                  } from '../subworkflows/local/delly'
include { DELLY_CNV              } from '../subworkflows/local/delly_cnv'
include { FETCH_REFERENCE_FILES  } from '../modules/local/fetch_reference_files'
include { COPY_BAMS              } from '../modules/local/copy_bams'
include { MERGE                  } from '../subworkflows/local/merge'

workflow CHORUS {
    take:
        ref_ch      // value channel: [ref_fa, ref_fai]
        fam_bam_ch  // queue: [fam, sam, bam, bai]

    main:
        vcfs = channel.empty()


        def needs_ref = params.callers.intersect(['SMOOVE', 'DELLY', 'DELLY_CNV'])
        if (needs_ref) { FETCH_REFERENCE_FILES() }
        
        if (params.copy_bams) { fam_bam_ch = COPY_BAMS(fam_bam_ch) }

        if (params.callers.contains('MANTA'))     { vcfs = vcfs.mix(MANTA(ref_ch, fam_bam_ch).vcfs) }
        if (params.callers.contains('SMOOVE'))    { vcfs = vcfs.mix(SMOOVE(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.smoove_excl).vcfs) }
        if (params.callers.contains('CNVNATOR'))  { vcfs = vcfs.mix(CNVNATOR(ref_ch, fam_bam_ch).vcfs) }
        if (params.callers.contains('DELLY'))     { vcfs = vcfs.mix(DELLY(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.delly_excl).vcfs) }
        if (params.callers.contains('DELLY_CNV')) { vcfs = vcfs.mix(DELLY_CNV(ref_ch, fam_bam_ch, FETCH_REFERENCE_FILES.out.delly_map).vcfs) }

        MERGE(vcfs)

    emit:
        vcfs      = vcfs                // [caller, sam, bcf, csi]
        collapsed = MERGE.out.collapsed // [sam, bcf, csi]
        merged    = MERGE.out.merged    // [cohort.bcf, cohort.bcf.csi]
}
