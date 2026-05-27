# nf-sv-pipe — CLAUDE.md

Nextflow SV-calling pipeline. Entry: `main.nf` → workflow `SVPLEX` in [workflows/sv_plex.nf](workflows/sv_plex.nf).

## Architecture

```
              ┌─ MATCHA  (per-sample collapse → [duphold] → cohort merge)   [params.matcha]
6 callers ───┼─ TRUVARI (per-sample collapse → [duphold] → cohort merge)   [params.truvari]
              └─ SVDB    (per-sample collapse → [duphold] → cohort merge)   [params.svdb]
```

Callers (MANTA, DYSGU, SMOOVE, CNVNATOR, DELLY, DELLY_CNV) each live in [subworkflows/local/](subworkflows/local/) and emit per-sample BCFs. All three merge branches consume the same post-`PASS_FILTER` channel and can be toggled independently via `params.matcha` / `params.truvari` / `params.svdb` (all default `true`). Merge logic is in [subworkflows/local/matcha.nf](subworkflows/local/matcha.nf), [subworkflows/local/truvari.nf](subworkflows/local/truvari.nf), and [subworkflows/local/svdb.nf](subworkflows/local/svdb.nf).

When `params.duphold` is true (default `true`), a [DUPHOLD](modules/local/duphold.nf) step runs between collapse and merge (applies to both branches). It restricts duphold annotation to DEL/DUP variants ≥ `params.duphold_min_size` (1 kb), then drops DELs with `FMT/DHFFC[0] > params.duphold_del_dhffc` (0.75) and DUPs with `FMT/DHBFC[0] < params.duphold_dup_dhbfc` (1.25). Smaller variants and other SVTYPE values bypass duphold entirely and are merged back before the cohort step. The `duphold` container label needs a combined bcftools+duphold image (see `nextflow.config`).

## Channel conventions

- Caller subworkflows emit `[caller, sam, bcf, csi]` — BCF + CSI, not VCF + TBI.
- `params.callers` list order encodes caller priority. `MERGE` re-sorts each per-sample tuple by `params.callers.indexOf(caller)` before handing to `MATCHA_COLLAPSE`; matcha keeps the first record on collapse.
- `ref_ch` is a value channel `[ref_fa, ref_fai]`. `fam_bam_ch` is a queue `[fam, sam, bam, bai]`.

## BCF/VCF conventions

- **Never pipe `bcftools concat` into `bcftools sort`**. To merge coordinate-sorted subsets use `bcftools concat --allow-overlaps` with tabix-indexed inputs.

## Tooling

- MATCHA modules use `label 'matcha'` — a combined matcha+bcftools container (`ghcr.io/jemunro/matcha/matcha-bcftools`, see `nextflow.config`). No separate binary is needed.
- TRUVARI modules use the `bcftools_truvari` label (a Seqera-built combined container with both tools). Per-sample collapse ([modules/local/truvari_collapse.nf](modules/local/truvari_collapse.nf)): strips all FORMAT fields except GT (`bcftools annotate -x '^FORMAT/GT'`) per-caller to avoid merge conflicts, then `bcftools merge -m id --force-samples` (same sample name across callers requires `--force-samples`; `-m id` is safe for symbolic alleles), splits into DEL/DUP/INV and BND/INS subsets, runs `truvari collapse --intra --chain` on each with type-specific params, then `bcftools concat | bcftools sort`. Cohort merge ([modules/local/truvari_merge.nf](modules/local/truvari_merge.nf)) follows the same split/typed-collapse pattern (without `--intra`) on the multi-sample merged VCF. The standalone `truvari` container label is also present but currently unused.
- DUPHOLD uses `label 'smoove'` — the smoove container bundles both `duphold` and `bcftools`.
- The `octopusv` and `jasminesv` container labels in [nextflow.config](nextflow.config) are dormant — no module currently references them. Same for the `jasmine_max_dist` param.

## Caller quirks

- **MANTA** emits a multi-sample family VCF; [modules/local/manta_split_sample.nf](modules/local/manta_split_sample.nf) splits to per-sample.
- **CNVNATOR** emits a text table; conversion to BCF goes via `cnvnator2VCF.awk` (in [bin/](bin/)) inside [modules/local/cnvnator_to_bcf.nf](modules/local/cnvnator_to_bcf.nf). Only `<DEL>` / `<DUP>` symbolic alleles.
- **DYSGU** calls per-sample only (no family-level merging). Emits a VCF; [modules/local/dysgu_to_bcf.nf](modules/local/dysgu_to_bcf.nf) converts to BCF.
- **DELLY_CNV** needs `params.delly_map` (mappability `.fa.gz`); it will fail without one.
- **SMOOVE / DELLY / DELLY_CNV** need reference auxiliary files, fetched on demand by `FETCH_REFERENCE_FILES` (gated on the caller list).

## Running tests

```
nextflow run . -profile test,singularity   # or test,conda
```

Uses the chr20 hg19 fixture in [test/fixtures/](test/fixtures/).
