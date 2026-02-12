#!/usr/bin/env Rscript

stopifnot(
    require(dplyr),
    require(matrixStats),
    require(stringr)
)

bin_fn   <- commandArgs(trailingOnly = TRUE)[1]
sm_cov   <- commandArgs(trailingOnly = TRUE)[-1]

bins <- readRDS(bin_fn)
samples <- sm_cov  %>% basename()  %>% str_remove('\\.coverage\\.rds$')
nsam <- length(sm_cov)
nbin <- nrow(bins)

COV <- matrix(0, nrow = nbin, ncol = nsam)
colnames(COV) <- samples

for (i in 1:nsam) {
    COV[,i] <- readRDS(sm_cov[i])
}

row_med <- rowMedians(COV, na.rm = TRUE)
row_mad <- rowMads(COV, constant = 1.4826, na.rm = TRUE)

CN <- 2 * sweep(COV, 1, row_med, "/")
Z <- sweep(COV, 1, row_med, "-") %>% sweep(1, row_mad, "/")

for (sm in samples) {
    bins %>%
        mutate(
            CN = CN[, sm],
            Z  = Z[, sm]
        ) %>%
        saveRDS(paste0(sm, '.bnorm.rds'))
}
