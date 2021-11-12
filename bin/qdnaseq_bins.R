#!/usr/bin/env Rscript
stopifnot(require(tidyverse),
          require(QDNAseq),
          require(docopt))

doc <- 
"
Usage:
  qdnaseq_bins.R <assembly> <output_bed>

Options:
  assembly        Assembly version, either hg19 or hg38
  output_bed      Name of output bed file
"

# opts <- docopt(doc, c('hg38', 'out.bed.gz'))
opts <- docopt(doc)

bins <- getBinAnnotations(binSize=10, genome = opts$assembly)
chr_pref <- `if`(opts$assembly == 'hg38', 'chr', '')

as_tibble(bins@data) %>% 
  select(chromosome, start, end) %>% 
  mutate(chromosome = str_c(chr_pref, chromosome)) %>% 
  write_tsv(opts$output_bed, col_names = F)

