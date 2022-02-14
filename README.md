# nf-cavalier

Nextflow cohort level SV calling pipeline based on Manta, Smoove, CNVnator and Jasmine. This pipeline is a work in progress.

## Usage
* Clone this repositoty
* Create and navigate to run directory
* See [bahlolab/nextflow-config](https://github.com/bahlolab/nextflow-config) for generic Nextflow configuration for Milton/SLURM.
* Create configuration file in run directory named `nextflow.config`, e.g.:
  ```Nextflow
    params {
      // inputs
      id = 'sv-run'
      ped = 'families.ped'
      bams = 'bams.tsv'
      
      // run config
      callers = ['MANTA', 'SMOOVE', 'CNVNATOR']
      assembly = 'hg38'
      ref_fasta = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/GATK/fasta_no_alt/hg38.no_alt.fasta'
    }
    ```
* **Params**  
  * `id` - Unique name for run. Used to name output files.
  * `ped` - Path to a [Ped format file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format) describing familial relationships, only first two columns used (family id and individual id).
  * `bams` - Path to TSV file with first column containing individual ID, second column containing path to indexed BAM file (no header row/  column names).
  * `callers` - List of SV callers to use. Currently supported values are "MANTA": [Manta](https://github.com/Illumina/manta), "CNVNATOR": [CNVnator](https://github.com/abyzovlab/CNVnator) and "SMOOVE": [Smoove](https://github.com/brentp/smoove) (lumpy).
* First run:  
`nextflow run /PATH/TO/nf-sv-pipe`
* Resume run:  
`nextflow run /PATH/TO/nf-sv-pipe -resume`
* Note: It is Recommended to run workflow in a `screen` session

## Output
* Outputs are created in the folder `output` in the run directory
* For each caller, a merged VCF file named `<id>.<caller>.vcf.gz` is output
* A combined VCF file including calls from all callers named `<id>.combined.vcf.gz` is also output

## Implementation
* Calling is implemented as recommended for individual callers
* Additional filtering for spilt-read/read-pair callers using [duphold](https://github.com/brentp/duphold) (Manta and Smoove only).
* Cohort merging using [Jasmine](https://github.com/mkirsche/Jasmine) (CNVnator and Manta only)
  * Breakend (SVTYPE=BND) calls are excluded at this point as they aren't handled well by Jasmine

## To-do
* Merging calls across different callers
* Implement additional callers (e.g. [GRIDDS](https://github.com/PapenfussLab/gridss))
