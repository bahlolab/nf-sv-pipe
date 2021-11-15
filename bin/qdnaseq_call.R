#!/usr/bin/env Rscript
stopifnot(require(tidyverse),
          require(QDNAseq),
          require(docopt))

doc <- 
"
Usage:
  qdnaseq_call.R <assembly> <ref_fai> <depth_bed> <name>

Options:
  assembly        Assembly version, either hg19 or hg38
  ref_fai         Reference fasta index for chromosome sizes
  depth_bed       BED file with depth column
  name            Name of sample, used for VCF
"

bin_read_counts <- function(bins, counts, sample_names=NULL) {
  
  if (is.null(sample_names)) {
    sample_names <- colnames(counts)
  }
  condition <- QDNAseq:::binsToUse(bins)
  phenodata <- 
    data.frame(name = sample_names,
               row.names = sample_names,
               paired.ends = TRUE,
               total.reads = colSums(counts),
               used.reads = colSums(counts[condition, , drop = FALSE]),
               stringsAsFactors = FALSE)
  
  object <- 
    new("QDNAseqReadCounts",
        bins = bins,
        counts = counts, 
        phenodata = phenodata)
  object$expected.variance <- QDNAseq:::expectedVariance(object)
  object
}

call_copy_num <- function(read_counts) {
  read_counts_filtered <- applyFilters(read_counts, residual=TRUE, blacklist=TRUE)
  read_counts_filtered <- estimateCorrection(read_counts_filtered)
  copy_numbers <- correctBins(read_counts_filtered)
  copy_numbers_normalized <- normalizeBins(copy_numbers)
  copy_numbers_smooth <- smoothOutlierBins(copy_numbers_normalized)
  copy_numbers_segmented <- segmentBins(copy_numbers_smooth)
  callBins(copy_numbers_segmented)
}

fix_vcf <- function(vcf_fn, bins, assembly, ref_fai) {
  
  chr_pref <- `if`(opts$assembly == 'hg38', 'chr', '')
  chrset <- str_c(chr_pref, unique(bins@data$chromosome))
  tmp_fn <- str_c(vcf_fn, '.tmp')
  
  contigs <-
    read_tsv(ref_fai,
             col_names = c('chrom', 'len', 'offset', 'lb', 'lw'),
             col_types = c('cidii')) %>% 
    filter(chrom %in% chrset) %>% 
    with(str_c('##contig=<ID=', chrom, ',length=', len,'>'))
  
  read_lines(vcf_fn) %>% 
    keep(~str_starts(., '##')) %>% 
    c("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">") %>% 
    c(contigs) %>% 
    write_lines(tmp_fn)
  
  read_tsv(vcf_fn, comment = '##', col_types = 'cccccccccc') %>% 
    mutate(`#CHROM` = str_c(chr_pref, `#CHROM`)) %>% 
    write_tsv(tmp_fn, append = TRUE, col_names = TRUE)
  
  file.rename(tmp_fn, vcf_fn)
  Rsamtools::bgzip(vcf_fn, str_c(vcf_fn, '.gz'))
}

# opts <- docopt(doc, c("hg38", "AH036.S35657_2.regions.bed.gz", "AH036.S35657_2"))
opts <- docopt(doc)

bins <- getBinAnnotations(binSize=10, genome = opts$assembly)

depth <- 
  read_tsv(opts$depth_bed, col_names = c('chrom', 'start', 'end', 'depth')) %>% 
  mutate(depth = as.integer(round(depth * 100)))

read_counts <- bin_read_counts(bins, matrix(depth$depth, ncol =1), sample_name = opts$name)
copy_num_calls <- call_copy_num(read_counts)
exportBins(copy_num_calls, format="vcf")

fix_vcf(str_c(opts$name, '.vcf'), bins, opts$assembly, opts$ref_fai)
