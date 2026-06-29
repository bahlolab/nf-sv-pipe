
include { read_caller_manifest  } from '../helpers'
include { read_merge_manifest   } from '../helpers'
include { build_caller_cache    } from '../helpers'
include { MANTA                 } from '../subworkflows/local/manta'
include { SMOOVE                } from '../subworkflows/local/smoove'
include { CNVNATOR              } from '../subworkflows/local/cnvnator'
include { DELLY                 } from '../subworkflows/local/delly'
include { DELLY_CNV             } from '../subworkflows/local/delly_cnv'
include { DYSGU                 } from '../subworkflows/local/dysgu'
include { FETCH_REFERENCE_FILES } from '../modules/local/fetch_reference_files'
include { MAKE_CALL_REGIONS     } from '../modules/local/make_call_regions'
include { COPY_BAMS             } from '../modules/local/copy_bams'
include { CRAM_TO_BAM           } from '../modules/local/cram_to_bam'
include { PASS_FILTER           } from '../modules/local/pass_filter'
include { MATCHA                } from '../subworkflows/local/matcha'
include { TRUVARI               } from '../subworkflows/local/truvari'
include { SVDB                  } from '../subworkflows/local/svdb'

workflow SVPLEX {
    take:
        ref_ch      // value channel: [ref_fa, [ref_fai, ref_gzi?]]
        chrs_ch     // value channel: List<String> (empty list = no restriction)
        fam_bam_ch  // queue: [fam, sam, bam, bai]; ignored in merge-only mode

    main:
        def outdir_abs = file(params.outdir).toAbsolutePath().toString()

        // ───── Merge-only mode ─────
        if (params.merge_manifest) {
            if (params.caller_manifest) log.warn "merge_manifest is set; ignoring params.caller_manifest"
            if (params.bams)            log.warn "merge_manifest is set; ignoring params.bams"
            if (params.ped)             log.warn "merge_manifest is set; ignoring params.ped"

            def merge_rows = read_merge_manifest(file(params.merge_manifest))
            def ch_for = { String br ->
                Channel.fromList(merge_rows.findAll { it.branch == br })
                       .map { r -> [r.sample, file(r.path), file(r.path + '.csi')] }
            }
            if (params.matcha)  MATCHA (channel.empty(), chrs_ch, channel.empty(), ref_ch, ch_for('MATCHA'))
            if (params.svdb)    SVDB   (channel.empty(), chrs_ch, channel.empty(), ref_ch, ch_for('SVDB'))
            if (params.truvari) TRUVARI(channel.empty(), chrs_ch, channel.empty(), ref_ch, ch_for('TRUVARI'))
            return
        }

        // ───── Normal mode ─────
        def caller_rows = params.caller_manifest ? read_caller_manifest(file(params.caller_manifest)) : []
        def cache       = build_caller_cache(caller_rows)
        def family_fully_cached = cache.family_fully_cached
        def sam_to_fam          = cache.sam_to_fam

        // Caller manifest entries usable as-is: family fully cached for the caller.
        def usable_caller_rows = caller_rows.findAll { r ->
            sam_to_fam.containsKey(r.sample) && sam_to_fam[r.sample] in family_fully_cached[r.caller]
        }

        fam_bam_ch = fam_bam_ch
            .groupTuple(by: 0)
            .map { fam, sam, bam, bai -> [groupKey(fam, sam.size()), sam, bam, bai] }
            .transpose()

        def needs_ref = params.callers.intersect(['MANTA', 'SMOOVE', 'DELLY', 'DELLY_CNV', 'CNVNATOR'])
        if (needs_ref) FETCH_REFERENCE_FILES()

        // Split CRAM from BAM; CRAMs always convert (implicit work-dir copy)
        split_input = fam_bam_ch.branch { _fam, _sam, bam, _bai ->
            cram: bam.name.endsWith('.cram')
            bam:  true
        }
        converted_ch = CRAM_TO_BAM(split_input.cram, ref_ch)

        // Optionally copy BAM inputs only (not CRAMs — already in work dir after conversion)
        def bam_ch
        if (params.copy_bams) {
            def will_run_any_caller = { fam -> params.callers.any { c -> !(fam.toString() in family_fully_cached[c]) } }
            def duphold_active      = params.duphold && (params.matcha || params.svdb || params.truvari)
            split_bams = split_input.bam.branch { fam, _sam, _bam, _bai ->
                to_copy:     duphold_active || will_run_any_caller(fam)
                passthrough: true
            }
            bam_ch = COPY_BAMS(split_bams.to_copy).mix(split_bams.passthrough)
        } else {
            bam_ch = split_input.bam
        }

        fam_bam_ch = bam_ch.mix(converted_ch)

        // Per-caller BAM input: drop families fully cached for that caller.
        def bam_for = params.callers.collectEntries { caller ->
            [(caller): fam_bam_ch.filter { fam, _sam, _bam, _bai -> !(fam.toString() in family_fully_cached[caller]) }]
        }

        vcfs_fresh = channel.empty()

        if (params.callers.contains('MANTA')) {
            call_regions_ch = MAKE_CALL_REGIONS(ref_ch, chrs_ch, FETCH_REFERENCE_FILES.out.delly_excl)
            MANTA(ref_ch, bam_for.MANTA, call_regions_ch)
            vcfs_fresh = vcfs_fresh.mix(MANTA.out.vcfs)
        }
        if (params.callers.contains('SMOOVE')) {
            SMOOVE(ref_ch, bam_for.SMOOVE, FETCH_REFERENCE_FILES.out.smoove_excl)
            vcfs_fresh = vcfs_fresh.mix(SMOOVE.out.vcfs)
        }
        if (params.callers.contains('CNVNATOR')) {
            CNVNATOR(ref_ch, chrs_ch, bam_for.CNVNATOR, FETCH_REFERENCE_FILES.out.delly_excl)
            vcfs_fresh = vcfs_fresh.mix(CNVNATOR.out.vcfs)
        }
        if (params.callers.contains('DELLY')) {
            DELLY(ref_ch, bam_for.DELLY, FETCH_REFERENCE_FILES.out.delly_excl)
            vcfs_fresh = vcfs_fresh.mix(DELLY.out.vcfs)
        }
        if (params.callers.contains('DELLY_CNV')) {
            DELLY_CNV(ref_ch, bam_for.DELLY_CNV, FETCH_REFERENCE_FILES.out.delly_map)
            vcfs_fresh = vcfs_fresh.mix(DELLY_CNV.out.vcfs)
        }
        if (params.callers.contains('DYSGU')) {
            DYSGU(ref_ch, bam_for.DYSGU)
            vcfs_fresh = vcfs_fresh.mix(DYSGU.out.vcfs)
        }

        manifest_caller_ch = Channel.fromList(usable_caller_rows)
            .map { r -> [r.caller, r.sample, file(r.path), file(r.path + '.csi')] }

        vcfs = vcfs_fresh.mix(manifest_caller_ch)

        branched = vcfs.branch { caller, _sam, _bcf, _csi ->
            filter:    params.apply_filters.contains(caller)
            no_filter: true
        }
        vcfs = PASS_FILTER(branched.filter).mix(branched.no_filter)

        sam_bam_ch = fam_bam_ch.map { _fam, sam, bam, bai -> [sam, bam, bai] }

        if (params.matcha)  MATCHA (vcfs, chrs_ch, sam_bam_ch, ref_ch, channel.empty())
        if (params.truvari) TRUVARI(vcfs, chrs_ch, sam_bam_ch, ref_ch, channel.empty())
        if (params.svdb)    SVDB   (vcfs, chrs_ch, sam_bam_ch, ref_ch, channel.empty())

        // ───── Output manifests (normal mode only) ─────
        // caller_manifest: fresh + usable passthrough (disjoint — fully-cached families never re-run).
        vcfs_fresh.map { caller, sam, bcf, _csi -> "${sam}\t${caller}\t${outdir_abs}/${caller}/${bcf.name}".toString() }
            .mix(Channel.fromList(usable_caller_rows).map { r -> "${r.sample}\t${r.caller}\t${r.path}".toString() })
            .collectFile(name: "${params.id}.caller_manifest.tsv", storeDir: params.outdir, newLine: true, sort: true)

        // merge_manifest: all fresh (merge-only mode does not pass through to this file).
        merge_rows_out = channel.empty()
        if (params.matcha)  merge_rows_out = merge_rows_out.mix(MATCHA.out.duphold .map { sam, bcf, _csi -> "${sam}\tMATCHA\t${outdir_abs}/MATCHA/${bcf.name}".toString()  })
        if (params.svdb)    merge_rows_out = merge_rows_out.mix(SVDB.out.duphold   .map { sam, bcf, _csi -> "${sam}\tSVDB\t${outdir_abs}/SVDB/${bcf.name}".toString()    })
        if (params.truvari) merge_rows_out = merge_rows_out.mix(TRUVARI.out.duphold.map { sam, bcf, _csi -> "${sam}\tTRUVARI\t${outdir_abs}/TRUVARI/${bcf.name}".toString() })
        merge_rows_out.collectFile(name: "${params.id}.merge_manifest.tsv", storeDir: params.outdir, newLine: true, sort: true)
}
