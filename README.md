# SVPLEX

Nextflow cohort-level SV calling pipeline using six callers ([MANTA](https://github.com/Illumina/manta), [DYSGU](https://github.com/kcleal/dysgu), [SMOOVE](https://github.com/brentp/smoove), [DELLY](https://github.com/dellytools/delly), [DELLY_CNV](https://github.com/dellytools/delly), [CNVNATOR](https://github.com/abyzovlab/CNVnator)) with cross-caller merging via [SVDB](https://github.com/J35P312/SVDB) (primary/default), [matcha](https://github.com/jemunro/matcha), and/or [truvari](https://github.com/acenglish/truvari).



## Usage

* Clone this repository
* Create and navigate to a run directory
* Create a `nextflow.config` in the run directory, e.g.:
  ```nextflow
  params {
    id        = 'sv-run'
    ped       = 'families.ped'
    bams      = 'bams.tsv'
    assembly  = 'hg38'
    ref_fasta = '/path/to/hg38.fasta'
    outdir    = 'output'
  }
  ```
* First run:
  ```
  nextflow run /PATH/TO/nf-sv-pipe
  ```
* Resume run:
  ```
  nextflow run /PATH/TO/nf-sv-pipe -resume
  ```
* Note: recommended to run in a `screen` or `tmux` session

## Params

| Param | Description |
|---|---|
| `id` | Unique name for the run; used in output filenames (default: `'SVPLEX'`) |
| `ped` | Path to a [PED format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format); only the first two columns (family ID, sample ID) are used |
| `bams` | Path to a TSV with sample ID in column 1 and path to indexed BAM/CRAM in column 2 (no header) |
| `ref_fasta` | Reference genome FASTA (must be indexed) |
| `assembly` | Genome build: `'hg38'` (default) or `'hg19'` |
| `outdir` | Output directory (default: `'output'`) |

<details>
<summary>Advanced parameters</summary>

> Note: the `delly_cnv_max_dels` / `delly_cnv_max_dups` and `duphold_max_dels` / `duphold_max_dups` defaults are empirically derived from outliers in a 30× 300+ paired-end WGS cohort — your mileage may vary; tune them for your data.

| Param | Description |
|---|---|
| `chrs` | Chromosomes to process; `null` = no restriction, `'auto'` = autosomes + X/Y (default), or a list of names |
| `callers` | List of callers to run; supported values: [`MANTA`](https://github.com/Illumina/manta), [`DYSGU`](https://github.com/kcleal/dysgu), [`SMOOVE`](https://github.com/brentp/smoove), [`DELLY`](https://github.com/dellytools/delly), `DELLY_CNV` (DELLY in CNV mode), [`CNVNATOR`](https://github.com/abyzovlab/CNVnator). Order sets merge priority (default: all six — `['MANTA', 'DYSGU', 'SMOOVE', 'DELLY', 'DELLY_CNV', 'CNVNATOR']`) |
| `svdb` | Run the SVDB merge branch — primary default (default: `true`) |
| `matcha` | Run the MATCHA merge branch (default: `false`) |
| `truvari` | Run the TRUVARI merge branch (default: `false`) |
| `apply_filters` | Callers whose BCFs are PASS-filtered before merging (default: `['DYSGU', 'DELLY']`) |
| `familial` | Group samples by family for joint calling where supported (default: `false`) |
| `duphold` | Run duphold between per-sample collapse and cohort merge to filter low-quality DEL/DUP calls (default: `true`) |
| `caller_manifest` | Optional TSV (`sample<TAB>caller<TAB>path`) of cached per-sample per-caller BCFs. When provided, caller calls are skipped for families whose every sample has a cached entry (joint-call all-or-nothing rule). Per-caller BCFs are auto-published to `${outdir}/<CALLER>/`. A new `${id}.caller_manifest.tsv` is written each run. Ignored in merge-only mode. |
| `merge_manifest` | Optional TSV (`sample<TAB>branch<TAB>path`) of per-sample post-collapse/duphold BCFs. **If set, triggers merge-only mode**: `params.bams` and `params.caller_manifest` are ignored (with a warning); only the samples listed are merged. Every sample must have an entry for every active merge branch. Per-branch DUPHOLD BCFs are auto-published to `${outdir}/<BRANCH>/`. A new `${id}.merge_manifest.tsv` is written in normal mode. |
| `chr_prefix` | Chromosome name prefix; `null` = auto-detect (`'chr'` for hg38, `''` for hg19) |
| `copy_bams` | Copy BAMs to work directory before calling — use when input is on slow or remote storage (default: `false`) |
| `refdir` | Directory for downloaded reference files (mappability, exclude lists); default: `'reference_files'` |
| `min_mapq` | Minimum mapping quality for reads (default: `15`) |
| `delly_cnv_max_dels` | If set, cap DELLY_CNV DEL callsets to top N by QUAL after CNV-normalisation (default: `2000`) |
| `delly_cnv_max_dups` | If set, cap DELLY_CNV DUP callsets to top N by QUAL after CNV-normalisation (default: `1000`) |
| `cnvnator_bin_size` | CNVnator bin size in bp (default: `1000`) |
| `cnvnator_exclude_overlap` | Drop CNVnator calls where the fraction of call length overlapping a delly exclude-list interval exceeds this threshold (default: `0.7`) |
| `duphold_min_size` | Minimum DEL/DUP size in bp for duphold annotation; smaller variants bypass duphold (default: `1000`) |
| `duphold_del_dhbfc` | Exclude DELs where `FMT/DHBFC[0]` exceeds this threshold (default: `0.75`) |
| `duphold_dup_dhbfc` | Exclude DUPs where `FMT/DHBFC[0]` is below this threshold (default: `1.25`) |
| `duphold_max_dels` | If set, cap DELs passing duphold by tightening DHBFC so at most N DELs pass (default: `4000`) |
| `duphold_max_dups` | If set, cap DUPs passing duphold by tightening DHBFC so at most N DUPs pass (default: `1000`) |
| `matcha_min_jaccard` | Minimum Jaccard similarity for matcha collapse/merge (default: `0.75`) |
| `matcha_sample_filter` | bcftools `-i` filter applied after per-sample matcha collapse (default: `'FILTER="PASS" || INFO/N_CALLERS>1'` — keeps PASS calls or events detected by ≥2 callers) |
| `matcha_cohort_filter` | bcftools `-i` filter applied after matcha cohort merge (default: `'INFO/N_CALLERS>1'` — keeps events detected by ≥2 callers in at least one sample) |
| `truvari_itvl_refdist` | DEL/DUP/INV collapse `--refdist` — max bp distance between breakpoints (default: `10000`) |
| `truvari_itvl_pctovl` | DEL/DUP/INV collapse `--pctovl` — min reciprocal overlap fraction (default: `0.75`) |
| `truvari_bnd_refdist` | BND/INS collapse `--refdist` — max bp distance between breakpoints (default: `50`) |
| `truvari_bnddist` | BND/INS collapse `--bnddist` — max bp distance for BND matching (default: `50`) |
| `truvari_bnd_pctsize` | BND/INS collapse `--pctsize` — min size similarity fraction (default: `0.75`) |
| `truvari_sample_filter` | bcftools `-i` filter applied after per-sample truvari collapse; default keeps events detected by ≥2 callers (`NumConsolidated>0`) |
| `truvari_cohort_filter` | bcftools `-i` filter applied after truvari cohort merge (default: `null` — no filter) |
| `svdb_overlap` | SVDB `--overlap`: min reciprocal overlap fraction for merging, both stages (default: `0.75`) |
| `svdb_bnd_distance` | SVDB `--bnd_distance`: max bp distance between precise breakpoints, both stages (default: `50`) |
| `svdb_sample_filter` | bcftools `-i` filter applied after per-sample SVDB collapse; default keeps events detected by ≥2 callers (`FOUNDBY>1`) |
| `svdb_cohort_filter` | bcftools `-i` filter applied after SVDB cohort merge (default: `null` — no filter) |
| `svdb_info_keep` | INFO fields retained in SVDB output BCFs (default: `['SVTYPE', 'SVLEN', 'END', 'POS2', 'CHR2', 'set', 'FOUNDBY', 'svdb_origin', 'SUPP_VEC']`) |

</details>

## Output

### Merged cohort BCFs (in `params.outdir`)

| File | Description |
|---|---|
| `<id>.MATCHA.cohort.bcf` (+ `.csi`) | Cohort-level BCF with per-sample matcha-collapsed calls merged across samples (only if `params.matcha`) |
| `<id>.TRUVARI.cohort.bcf` (+ `.csi`) | Cohort-level BCF with per-sample truvari-collapsed calls merged across samples (only if `params.truvari`) |
| `<id>.SVDB.cohort.bcf` (+ `.csi`) | Cohort-level BCF with per-sample SVDB-collapsed calls merged across samples (only if `params.svdb`) |

### Per-sample collapsed BCFs (in `params.outdir`)

| File | Description |
|---|---|
| `MATCHA/<sample>.MATCHA.bcf` (+ `.csi`) | Per-sample BCF with calls from all callers collapsed by matcha (only if `params.matcha`) |
| `TRUVARI/<sample>.TRUVARI.bcf` (+ `.csi`) | Per-sample BCF with calls from all callers collapsed by truvari (only if `params.truvari`) |
| `SVDB/<sample>.SVDB.bcf` (+ `.csi`) | Per-sample BCF with calls from all callers collapsed by SVDB (only if `params.svdb`) |

### Per-caller per-sample BCFs (in `${outdir}/<CALLER>/`)

The raw per-sample BCF from each caller is auto-published under `${params.outdir}/<CALLER>/`. Filename convention: `<sample>.<CALLER>.bcf` (+ `.csi`). For DELLY samples in multi-member families the file is `<sample>.DELLY.geno.bcf` (genotyped against a joint site list). These paths are listed in `${id}.caller_manifest.tsv` and can be fed back into a later run via `--caller_manifest`.

### Cross-batch merge inputs (in `${outdir}/<BRANCH>/`)

The per-sample post-collapse/duphold BCF for each active merge branch (MATCHA/SVDB/TRUVARI) is auto-published under `${params.outdir}/<BRANCH>/`. These paths are listed in `${id}.merge_manifest.tsv`; concatenating per-batch merge manifests and passing the result via `--merge_manifest` triggers merge-only mode (no calls run, only cohort merging).

## Manifests

Every normal run writes two manifest TSVs alongside the cohort outputs:

| File | Columns | What it lists | Use it as input via |
|---|---|---|---|
| `${id}.caller_manifest.tsv` | `sample`, `caller`, `path` | Per-sample raw caller BCFs published under `${outdir}/<CALLER>/` | `--caller_manifest` |
| `${id}.merge_manifest.tsv`  | `sample`, `branch`, `path` | Per-sample post-collapse/duphold BCFs (one per active merge branch) published under `${outdir}/<BRANCH>/` | `--merge_manifest` |

The two manifests serve **independent** use cases:

### `caller_manifest` — cache calling

When a `caller_manifest` is provided, the pipeline skips re-calling for any sample with a cached entry for a given caller. Everything downstream (PASS filter, collapse, duphold, cohort merge) still runs.

Typical use: persistent caching of caller calls. Nextflow's own `-resume` cache lives in the work directory, which on shared HPC scratch is often auto-cleaned; the `caller_manifest` lets you skip re-running callers on a fresh work dir.

```bash
# First run — produces out/run1.caller_manifest.tsv
nextflow run /PATH/TO/nf-sv-pipe --bams cohort.bams --ped fam.ped \
    --outdir out --id run1

# Later run (e.g. after scratch cleanup) — reuses cached caller BCFs
nextflow run /PATH/TO/nf-sv-pipe --bams cohort.bams --ped fam.ped \
    --caller_manifest out/run1.caller_manifest.tsv \
    --outdir out --id run2
```

(When `params.familial=true`, joint callers MANTA/DELLY are skipped per family only if every sample in the family has a cached entry; otherwise the caller re-runs for the whole family.)

### `merge_manifest` — merge across batches processed separately

When all BAMs can't fit on disk simultaneously, split the cohort into batches, run each batch independently (each writes its own `merge_manifest`), then run **once more in merge-only mode** with the concatenated manifests to produce the cohort:

```bash
# Batch 1
nextflow run /PATH/TO/nf-sv-pipe --bams batch1.bams --ped fam.ped \
    --outdir out_b1 --id batch1
# → out_b1/batch1.merge_manifest.tsv

# Batch 2
nextflow run /PATH/TO/nf-sv-pipe --bams batch2.bams --ped fam.ped \
    --outdir out_b2 --id batch2

# Concatenate the per-batch merge manifests
cat out_b{1,2}/*.merge_manifest.tsv > combined.merge_manifest.tsv

# Cross-batch merge — runs ONLY cohort merging. No BAMs needed.
nextflow run /PATH/TO/nf-sv-pipe \
    --merge_manifest combined.merge_manifest.tsv \
    --outdir out_merge --id merged
```

Setting `--merge_manifest` puts the pipeline into **merge-only mode**: callers, PASS filter, per-sample collapse, and duphold are *all* skipped. Any `--bams` / `--ped` / `--caller_manifest` supplied alongside is ignored (with a warning). Every sample in the manifest must have an entry for every active merge branch (MATCHA / SVDB / TRUVARI per `params.matcha`/`svdb`/`truvari`).

> **Note**: the multi-batch workflow requires `params.duphold=true` (the default). With `--duphold false`, `merge_manifest.tsv` is written empty because the per-sample BCFs fed to cohort merge are COLLAPSE outputs, which are not published to a durable path. Cross-batch merging in that configuration is not currently supported.

## Implementation

* Each enabled caller runs in parallel. Callers that support joint family calling (MANTA, SMOOVE, DELLY, DELLY_CNV) operate on family-grouped BAMs.
* Callers listed in `apply_filters` (default: DYSGU, DELLY) have their per-sample BCFs PASS-filtered before merging.
* The PASS-filtered caller channel feeds all three merge branches; each runs independently and can be toggled with `params.svdb`, `params.matcha`, or `params.truvari`. SVDB is the primary default (`true`); MATCHA and TRUVARI default to `false`.
* **MATCHA per-sample collapse**: `matcha collapse` merges calls from all callers per sample. Caller priority follows `params.callers` order — the first caller's record is kept when calls are collapsed.
* **TRUVARI per-sample collapse**: FORMAT fields are stripped to GT-only per caller, then `bcftools merge -m id --force-samples` creates a multi-column VCF (one column per caller). Variants are split into DEL/DUP/INV (interval-overlap matching) and BND/INS (breakpoint-proximity matching) subsets; `truvari collapse --intra --chain` runs on each with type-specific params, and results are concatenated and sorted. `--intra` consolidates the per-caller columns into a single sample column with a `FORMAT/SUPP` support field.
* **DUPHOLD** (when `params.duphold = true`): After per-sample collapse, duphold annotates DEL/DUP variants ≥ `params.duphold_min_size` with depth-fold-change fields (`FMT/DHBFC`). DELs with `DHBFC > params.duphold_del_dhbfc` and DUPs with `DHBFC < params.duphold_dup_dhbfc` are dropped. If `params.duphold_max_dels` / `params.duphold_max_dups` are set and more than N DELs / DUPs would pass, the DHBFC threshold is tightened dynamically — the Nth-best metric is used in place of the default, so at most N records of that type survive. Smaller variants and non-DEL/DUP SVTYPE values bypass duphold and are merged back before the cohort step. Applies to all merge branches.
* **DELLY_CNV normalisation**: DELLY_CNV emits `SVTYPE=CNV` records, which are recoded to DEL/DUP based on copy number by `bin/delly_cnv_norm.awk`. If `params.delly_cnv_max_dels` / `params.delly_cnv_max_dups` are set, the normalised callset is then capped per-type to the top N by QUAL.
* **MATCHA cohort merge**: `matcha merge` pools all per-sample collapsed BCFs into a single cohort BCF.
* **TRUVARI cohort merge**: `bcftools merge -m id` pools per-sample collapsed BCFs into a multi-sample VCF, then the same DEL/DUP/INV / BND/INS split-collapse approach runs without `--intra`.
* **SVDB per-sample collapse**: BCFs are converted to VCF.gz; `svdb --merge --priority <callers>` merges them with caller priority matching `params.callers` order.
* **SVDB cohort merge**: Per-sample SVDB-collapsed BCFs are converted to VCF.gz; `svdb --merge` merges across samples. No caller priority is applied at this stage.
