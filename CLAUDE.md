# nf-sv-pipe — CLAUDE.md

Nextflow SV-calling pipeline. Entry: `main.nf` → workflow `SVPLEX` in [workflows/sv_plex.nf](workflows/sv_plex.nf).

## Architecture

```
6 callers → MATCHA_COLLAPSE (per-sample) → MATCHA_MERGE (cohort)
```

Callers (MANTA, DYSGU, SMOOVE, CNVNATOR, DELLY, DELLY_CNV) each live in [subworkflows/local/](subworkflows/local/) and emit per-sample BCFs. Merge logic is in [subworkflows/local/matcha.nf](subworkflows/local/matcha.nf).

## Channel conventions

- Caller subworkflows emit `[caller, sam, bcf, csi]` — BCF + CSI, not VCF + TBI.
- `params.callers` list order encodes caller priority. `MERGE` re-sorts each per-sample tuple by `params.callers.indexOf(caller)` before handing to `MATCHA_COLLAPSE`; matcha keeps the first record on collapse.
- `ref_ch` is a value channel `[ref_fa, ref_fai]`. `fam_bam_ch` is a queue `[fam, sam, bam, bai]`.

## Tooling

- `matcha` binary is expected at `bin/matcha` but is **not committed to git** — it must be obtained/built separately. Its modules use `container null` and invoke it via `bin/` (Nextflow adds `bin/` to PATH automatically for processes in this repo).
- The `octopusv`, `truvari`, and `jasminesv` container labels in [nextflow.config](nextflow.config) are dormant — no module currently references them. Same for the `truvari_*` and `jasmine_max_dist` params.

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
