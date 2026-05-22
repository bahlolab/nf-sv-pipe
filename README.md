# nf-sv-plex

Nextflow cohort-level SV calling pipeline using six callers ([MANTA](https://github.com/Illumina/manta), [DYSGU](https://github.com/kcleal/dysgu), [SMOOVE](https://github.com/brentp/smoove), [DELLY](https://github.com/dellytools/delly), DELLY_CNV, [CNVNATOR](https://github.com/abyzovlab/CNVnator)) with cross-caller merging via matcha and/or [truvari](https://github.com/acenglish/truvari) (toggle either branch with `params.matcha` / `params.truvari`).

## Prerequisites

* `matcha` binary must be placed at `bin/matcha` in the pipeline directory before running. It is not distributed with this repo.
* See [bahlolab/nextflow-config](https://github.com/bahlolab/nextflow-config) for generic Nextflow configuration for Milton/SLURM.

## Usage

* Clone this repository
* Create and navigate to a run directory
* Create a `nextflow.config` in the run directory, e.g.:
  ```nextflow
  params {
    // inputs
    id        = 'sv-run'
    ped       = 'families.ped'
    bams      = 'bams.tsv'

    // run config
    callers   = ['MANTA', 'DYSGU', 'SMOOVE', 'DELLY', 'DELLY_CNV', 'CNVNATOR']
    assembly  = 'hg38'
    ref_fasta = '/path/to/hg38.fasta'
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
| `id` | Unique name for the run; used in output filenames |
| `ped` | Path to a [PED format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format); only the first two columns (family ID, sample ID) are used |
| `bams` | Path to a TSV with sample ID in column 1 and path to indexed BAM/CRAM in column 2 (no header) |
| `ref_fasta` | Reference genome FASTA (must be indexed) |
| `assembly` | Genome build: `'hg38'` (default) or `'hg19'` |
| `callers` | List of callers to run; supported values: [`MANTA`](https://github.com/Illumina/manta), [`DYSGU`](https://github.com/kcleal/dysgu), [`SMOOVE`](https://github.com/brentp/smoove), [`DELLY`](https://github.com/dellytools/delly), `DELLY_CNV` (DELLY in CNV mode), [`CNVNATOR`](https://github.com/abyzovlab/CNVnator). Order sets merge priority. |
| `apply_filters` | Callers whose BCFs are PASS-filtered before merging (default: `['DYSGU', 'DELLY']`) |
| `familial` | Group samples by family for joint calling where supported (default: `true`) |
| `chr_prefix` | Chromosome name prefix; `null` = auto-detect (`'chr'` for hg38, `''` for hg19) |
| `chrs` | Chromosomes to process; `null` = no restriction, `'auto'` = autosomes + X/Y (default), or a list of names |
| `copy_bams` | Copy BAMs to work directory before calling — use when input is on slow or remote storage (default: `false`) |
| `refdir` | Directory for downloaded reference files (mappability, exclude lists); default: `'reference_files'` |
| `cachedir` | `storeDir` path for cacheable call outputs; `null` = always re-run (default) |
| `cnvnator_bin_size` | CNVnator bin size in bp (default: `1000`) |
| `matcha` | Run the MATCHA merge branch (default: `true`) |
| `matcha_min_jaccard` | Minimum Jaccard similarity for matcha collapse/merge (default: `0.75`) |
| `matcha_sample_filter` | bcftools filter expression applied after per-sample matcha collapse; default keeps PASS or multi-caller calls |
| `matcha_cohort_filter` | bcftools filter expression applied after matcha cohort merge; default keeps multi-caller calls |
| `truvari` | Run the TRUVARI merge branch (default: `true`) |
| `truvari_itvl_refdist` | DEL/DUP/INV collapse `--refdist` — max bp distance between breakpoints (default: `10000`) |
| `truvari_itvl_pctovl` | DEL/DUP/INV collapse `--pctovl` — min reciprocal overlap fraction (default: `0.75`) |
| `truvari_bnd_refdist` | BND/INS collapse `--refdist` — max bp distance between breakpoints (default: `50`) |
| `truvari_bnddist` | BND/INS collapse `--bnddist` — max bp distance for BND matching (default: `50`) |
| `truvari_bnd_pctsize` | BND/INS collapse `--pctsize` — min size similarity fraction (default: `0.75`) |
| `truvari_sample_filter` | bcftools filter expression applied after per-sample truvari collapse; default keeps PASS or consolidated calls |
| `truvari_cohort_filter` | bcftools filter expression applied after truvari cohort merge (default: `null` — no filter) |

## Output

Outputs are written to `params.outdir` (default: `output/`) in the run directory:

| File | Description |
|---|---|
| `collapse/<sample>.collapsed.bcf` (+ `.csi`) | Per-sample BCF with calls from all callers collapsed by matcha (only if `params.matcha`) |
| `<id>.cohort.bcf` (+ `.csi`) | Cohort-level BCF with per-sample matcha-collapsed calls merged across samples (only if `params.matcha`) |
| `truvari_collapse/<sample>.truvari.collapsed.bcf` (+ `.csi`) | Per-sample BCF with calls from all callers collapsed by truvari (only if `params.truvari`) |
| `<id>.truvari.cohort.bcf` (+ `.csi`) | Cohort-level BCF with per-sample truvari-collapsed calls merged across samples (only if `params.truvari`) |

## Implementation

* Each enabled caller runs in parallel. Callers that support joint family calling (MANTA, SMOOVE, DELLY, DELLY_CNV) operate on family-grouped BAMs.
* Callers listed in `apply_filters` (default: DYSGU, DELLY) have their per-sample BCFs PASS-filtered before merging.
* The PASS-filtered caller channel feeds **both** the MATCHA and TRUVARI branches; each branch runs independently and can be disabled with `params.matcha = false` or `params.truvari = false`.
* **MATCHA per-sample collapse**: `matcha collapse` merges calls from all callers per sample. Caller priority follows `params.callers` order — the first caller's record is kept when calls are collapsed.
* **MATCHA cohort merge**: `matcha merge` pools all per-sample collapsed BCFs into a single cohort BCF.
* **TRUVARI per-sample collapse**: FORMAT fields are stripped to GT-only per caller, then `bcftools merge -m id --force-samples` creates a multi-column VCF (one column per caller). Variants are split into DEL/DUP/INV (interval-overlap matching) and BND/INS (breakpoint-proximity matching) subsets; `truvari collapse --intra --chain` runs on each with type-specific params, and results are concatenated and sorted. `--intra` consolidates the per-caller columns into a single sample column with a `FORMAT/SUPP` support field.
* **TRUVARI cohort merge**: `bcftools merge -m id` pools per-sample collapsed BCFs into a multi-sample VCF, then the same DEL/DUP/INV / BND/INS split-collapse approach runs without `--intra`.
