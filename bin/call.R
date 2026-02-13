#!/usr/bin/env Rscript

stopifnot(
    require(readr),
    require(dplyr),
    require(tidyr),
    require(stringr),
    require(purrr)
)

output <- commandArgs(trailingOnly = TRUE)[1]
inputs <- commandArgs(trailingOnly = TRUE)[-1]

MIN_RUN <- 5L
MAX_GAP <- 1L
MIN_P   <- 0.01
MAX_DEL <- 1.5
MIN_DUP <- 2.5
MAX_DUP <- 6

MAX_GAP_FRAC <- 0.25
#

bnorm <- map_df(inputs, readRDS)

# call per phase
calls_1 <-
    bnorm %>%
    filter(
        !is.na(CN),
    ) %>% 
    mutate(
        IDX = row_number(),
        B_PVAL = 2 * pnorm(-abs(Z))
    ) %>%
    mutate(
        ZZ = (Z - median(Z, na.rm = TRUE)) / mad(Z, na.rm = TRUE),
        S_PVAL = 2 * pnorm(-abs(ZZ))
    ) %>%
    ungroup() %>%
    filter(
        B_PVAL < MIN_P,
        S_PVAL < MIN_P,
        CN < MAX_DEL | (CN > MIN_DUP & CN < MAX_DUP)
    ) %>%
    group_by(chrom, phase) %>% 
    mutate(
        STEP = IDX - lag(IDX),
        BREAK = is.na(STEP) | STEP <= 0 | STEP > 1 + MAX_GAP,
        SEG = cumsum(BREAK)
    ) %>%
    group_by(chrom, phase, SEG) %>%
    filter(n() >= MIN_RUN)  %>% 
    summarise(
        start = first(start),
        end = last(end),
        CN = median(CN),
        SVTYPE = if_else(CN > 2, "DUP", "DEL"),
        .groups = "drop"
    ) %>%
    select(-SEG)

# merge phases
calls_2 <-
    calls_1 %>%
    arrange(SVTYPE, chrom, start, end) %>%
    group_by(SVTYPE, chrom) %>%
    mutate(
        BREAK = is.na(lag(end)) | start > lag(end),
        SEG = cumsum(BREAK)
    ) %>%
    group_by(SVTYPE, chrom, SEG, phase) %>%
    summarise(
        CN = weighted.mean(CN, end - start),
        start = min(start),
        end = max(end),
        .groups = "drop"
    ) %>%
    group_by(SVTYPE, chrom, SEG) %>%
    summarise(
        start = median(start),
        end = median(end),
        CN = median(CN),
        SVLEN = end - start,
        .groups = "drop"
    ) %>%
    select(-SEG)

# merge across gaps
calls_3 <-
    calls_2 %>%
    arrange(SVTYPE, chrom, start) %>%
    group_by(SVTYPE, chrom) %>%
    (function(calls) {
        while (TRUE) {
            new_calls <-
                calls %>%
                mutate(
                    GAP_FRAC = (start - lag(end)) / (end - lag(start)),
                    BREAK = is.na(lag(end)) | GAP_FRAC > MAX_GAP_FRAC,
                    SEG = cumsum(BREAK)
                ) %>%
                group_by(SVTYPE, chrom, SEG) %>%
                summarise(
                    start = min(start),
                    end = max(end),
                    CN = weighted.mean(CN, SVLEN, na.rm = TRUE),
                    SVLEN = end - start,
                    .groups = "drop"
                ) %>%
                select(-SEG)

            if (nrow(calls) == nrow(new_calls)) {
                return(calls)
            }
            calls <- new_calls
        }
    }) %>%
    ungroup() %>%
    mutate(
        GT = case_when(
            CN < 0.5 ~ "1/1",
            CN < 3.5 ~ "0/1",
            TRUE ~ "1/1",
        )
    )  %>% 
    arrange(chrom, start, end)

vcf_lines <- c(
  "##fileformat=VCFv4.2",
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
  '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number estimate">',
  str_c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", output),
  calls_3 %>%
    transmute(
      CHROM  = chrom,
      POS    = as.integer(start) + 1L,
      ID     = '.',
      REF    = "N",
      ALT    = str_c("<", SVTYPE, ">"),
      QUAL   = ".",
      FILTER = "PASS",
      INFO   = str_c("SVTYPE=", SVTYPE, ";END=", end, ";SVLEN=", if_else(SVTYPE == "DEL", -SVLEN, SVLEN)),
      FORMAT = "GT:CN",
      SAMPLE = str_c(GT, ":", formatC(CN, format = "f", digits = 3))
    ) %>%
    mutate(line = str_c(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE, sep = "\t")) %>%
    pull(line)
)

write_lines(vcf_lines, str_c(output, '.vcf'))
