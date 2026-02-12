#!/usr/bin/env Rscript

stopifnot(
    require(readr),
    require(dplyr)
)

input  <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
MIN_RUN <- 5L
MAX_GAP <- 1L
MIN_P   <- 0.05
MAX_DEL <- 1.5
MIN_DUP <- 2.5
MAX_DUP <- 6

bnorm <- readRDS(input)

bnorm %>%
    mutate(
        P = 2 * pnorm(-abs(Z)),
        IDX = row_number()
    ) %>%
    filter(
        !is.na(CN),
        P < MIN_P,
        CN < MAX_DEL | (CN > MIN_DUP & CN < MAX_DUP)
    ) %>%
    mutate(
        STEP = IDX - lag(IDX),
        STEP = if_else(chr != lag(chr), Inf, STEP),
        BREAK = is.na(STEP) | STEP <= 0 | STEP > 1 + MAX_GAP,
        SEG = cumsum(BREAK)
    ) %>%
    group_by(SEG) %>%
    filter(n() >= MIN_RUN) %>%
    summarise(
        CHROM = first(chr),
        POS = first(start),
        END = last(end),
        SVLEN = END - POS,
        CN = median(CN),
        SVTYPE = if_else(CN > 2, "DUP", "DEL"),
        GT = case_when(
            CN < 0.5 ~ "1/1",
            CN < 3.5 ~ "0/1",
            TRUE ~ "1/1",
        ),
        .groups = "drop"
    ) %>%
    select(-SEG) %>%
    write_tsv(paste0(output, ".CNVs.tsv"))
