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

# modified from QDNAseq::exportVCF to include FC fold change annotation and filtering
export_vcf <- function(obj, counts_obj) {
  counts <- Biobase::assayDataElement(counts_obj, "counts")
  calls <- Biobase::assayDataElement(obj, "calls")
  segments <- QDNAseq:::log2adhoc(Biobase::assayDataElement(obj, "segmented"))
  fd <- Biobase::fData(obj)
  pd <- Biobase::pData(obj)
  vcfHeader <- cbind(c("##fileformat=VCFv4.2", paste("##source=QDNAseq-", 
                                                     packageVersion("QDNAseq"), sep = ""), "##REF=<ID=DIP,Description=\"CNV call\">", 
                       "##ALT=<ID=DEL,Description=\"Deletion\">", "##ALT=<ID=DUP,Description=\"Duplication\">", 
                       "##FILTER=<ID=DEL_HIGH_FC,Description=\"Filtered due to high fold change deletion\">", 
                       "##FILTER=<ID=DUP_LOW_FC,Description=\"Filtered due to low fold change duplication\">", 
                       "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant: DEL,DUP,INS\">", 
                       "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of variant\">", 
                       "##INFO=<ID=BINS,Number=1,Type=Integer,Description=\"Number of bins in call\">", 
                       "##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Score of calling algorithm\">", 
                       "##INFO=<ID=LOG2CNT,Number=1,Type=Float,Description=\"Log 2 count\">", 
                       "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                       "##FORMAT=<ID=FC,Number=1,Type=Float,Description=\"Raw fold change vs median depth\">"))
  for (i in 1:ncol(calls)) {
    d <- cbind(fd[, 1:3], calls[, i], segments[, i])
    sel <- d[, 4] != 0 & !is.na(d[, 4])
    ### calc median depth
    inv_sel <- d[, 4] == 0 & !is.na(d[, 4])
    median_depth <- median(counts[inv_sel, 1])
    ###
    dsel <- d[sel, ]
    counts_sel <- counts[sel, i]
    rleD <- rle(paste(d[sel, 1], d[sel, 4], sep = ":"))
    endI <- cumsum(rleD$lengths)
    posI <- c(1, endI[-length(endI)] + 1)
    fc <- map2_dbl(posI, endI, function(p, e) { mean(counts_sel[p:e]) / median_depth })
    chr <- dsel[posI, 1]
    pos <- dsel[posI, 2]
    end <- dsel[endI, 3]
    score <- dsel[posI, 4]
    segVal <- round(dsel[posI, 5], 2)
    svtype <- rep(NA, length(chr))
    svlen <- rep(NA, length(chr))
    gt <- rep(NA, length(chr))
    bins <- rleD$lengths
    svtype[dsel[posI, 4] <= -1] <- "DEL"
    svtype[dsel[posI, 4] >= 1] <- "DUP"
    svlen <- end - pos + 1
    gt[score == -2] <- "1/1"
    gt[score == -1] <- "0/1"
    gt[score == 1] <- "0/1"
    gt[score == 2] <- "0/1"
    gt[score == 3] <- "0/1"
    options(scipen = 100)
    id <- "."
    ref <- "<DIP>"
    alt <- paste("<", svtype, ">", sep = "")
    qual <- 1000
    filter <- "PASS"
    filter <- case_when(svtype == 'DEL' & fc > 0.75 ~ 'DEL_HIGH_FC',
                        svtype == 'DUP' & fc < 1.25 ~ 'DUP_LOW_FC',
                        TRUE ~ 'PASS')
    info <- paste("SVTYPE=", svtype, ";END=", end, ";SVLEN=", 
                  svlen, ";BINS=", bins, ";SCORE=", score, ";LOG2CNT=", 
                  segVal, sep = "")
    format <- "GT:FC"
    sample <- str_c(gt, ':', format(round(fc, 3)))
    out <- cbind(chr, pos, id, ref, alt, qual, filter, info, 
                 format, sample)
    colnames(out) <- c("#CHROM", "POS", "ID", "REF", "ALT", 
                       "QUAL", "FILTER", "INFO", "FORMAT", pd$name[i])
    fname <- paste(pd$name[i], ".vcf", sep = "")
    write.table(vcfHeader, fname, quote = F, sep = "\t", 
                col.names = FALSE, row.names = FALSE)
    suppressWarnings(write.table(out, fname, quote = F, 
                                 sep = "\t", append = TRUE, col.names = TRUE, row.names = FALSE))
  }
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

# opts <- docopt(doc, c('hg38', 'hg38.no_alt.fasta.fai', 'ST734_3.regions.bed.gz', 'ST734_3'))
opts <- docopt(doc)

bins <- getBinAnnotations(binSize=10, genome = opts$assembly)

depth <- 
  read_tsv(opts$depth_bed, col_names = c('chrom', 'start', 'end', 'depth')) %>% 
  mutate(depth = as.integer(round(depth * 100)))

read_counts <- bin_read_counts(bins, matrix(depth$depth, ncol =1), sample_name = opts$name)

copy_num_calls <- call_copy_num(read_counts)

saveRDS(copy_num_calls, str_c(opts$name, '.rds'))

export_vcf(copy_num_calls, read_counts)

fix_vcf(str_c(opts$name, '.vcf'), bins, opts$assembly, opts$ref_fai)
