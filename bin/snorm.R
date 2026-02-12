#!/usr/bin/env Rscript

stopifnot(
    require(readr),
    require(dplyr)
)

input  <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]

cov <-
    read_tsv(
        input,
        col_names = c("chr", "start", "end", "coverage"),
        col_types = "ciin"
    ) %>%
    filter(chr %in% paste0("chr", c(1:22))) %>%
    mutate(
        coverage = coverage / median(
            coverage[coverage > 0],
            na.rm = TRUE
        )
    )

saveRDS(cov$coverage, paste0(output, '.coverage.rds'))

cov %>%
    select(-coverage) %>%
    saveRDS(paste0(output, '.bins.rds'))

