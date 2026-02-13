#!/usr/bin/env Rscript

stopifnot(
    require(dplyr),
    require(matrixStats),
    require(stringr)
)

bin_fn  <- commandArgs(trailingOnly = TRUE)[1]
snorm   <- commandArgs(trailingOnly = TRUE)[-1]

bins <- readRDS(bin_fn)
samples <- snorm  %>% basename() %>% str_remove('\\.snorm\\.rds$')
nsam <- length(snorm)
nbin <- nrow(bins)

SNORM <- matrix(0, nrow = nbin, ncol = nsam)
colnames(SNORM) <- samples

for (i in 1:nsam) {
    SNORM[,i] <- readRDS(snorm[i])
}

row_med <- rowMedians(SNORM, na.rm = TRUE)
row_mad <- rowMads(SNORM, constant = 1.4826, na.rm = TRUE)

CN <- 2 * sweep(SNORM, 1, row_med, "/")
Z <- sweep(SNORM, 1, row_med, "-") %>% sweep(1, row_mad, "/")

for (sm in samples) {
    bins %>%
        mutate(
            CN = CN[, sm],
            Z  = Z[, sm]
        ) %>%
        saveRDS(paste0(sm, '.bnorm.rds'))
}
