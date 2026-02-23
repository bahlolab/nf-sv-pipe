#!/usr/bin/env Rscript

stopifnot(
    require(readr),
    require(dplyr),
    require(purrr)
)

input  <- commandArgs(trailingOnly = TRUE)[1]
output <- commandArgs(trailingOnly = TRUE)[2]
NPHASE <- as.integer(commandArgs(trailingOnly = TRUE)[3])
NSHARDS <- as.integer(commandArgs(trailingOnly = TRUE)[4])
NUC <- commandArgs(trailingOnly = TRUE)[5]

MAX_PN <- 0.25

pool_phase_mean <- function(df, NPHASE = 2L) {

  NPHASE <- as.integer(NPHASE)

  df %>%
    arrange(chrom, start, end) %>%
    group_by(chrom) %>%
    reframe({
      n0 <- n()
      if (n0 < NPHASE) return(tibble())

      v  <- .data$depth
      cs <- c(0, cumsum(v))

      bind_rows(lapply(seq_len(NPHASE) - 1L, function(off) {

        starts <- seq.int(1L + off, n0 - NPHASE + 1L, by = NPHASE)
        if (!length(starts)) return(NULL)

        ends <- starts + NPHASE - 1L

        tibble(
          start = start[starts],
          end   = end[ends],
          depth = (cs[ends + 1L] - cs[starts]) / NPHASE,
          phase = off + 1L
        )
      }))
    }) %>%
    ungroup()
}

gc_depth <-
  read_tsv(NUC) %>%
  filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
  mutate(
    nACTG = nA + nC + nG + nT,
    pN = nN / (nN + nACTG),
    pGC = (nG + nC) / nACTG
  ) %>%
  filter(pN < MAX_PN) %>%
  select(chrom, start, end, pGC) %>%
  inner_join(
    read_tsv(
      input,
      col_names = c("chrom", "start", "end", "depth"),
      col_types = "ciin"
    ),
    by = join_by(chrom, start, end)
  ) %>% 
  group_by(chrgrp = chrom == "chrX") %>%
  mutate(
      depth = depth / median(
        depth[depth > 0],
        na.rm = TRUE
      )
  ) %>%
  ungroup()  %>% 
  select(-chrgrp)

gc_binned <-
  gc_depth %>%
  filter(depth > 0) %>%
  mutate(
    ld = log(depth),
    z = (ld - median(ld)) / mad(ld)
  ) %>%
  filter(abs(z) < -qnorm(0.05)) %>%
  arrange(pGC) %>%
  mutate(bin = cut(row_number(), 10000, labels = seq_len(10000))) %>%
  group_by(bin) %>%
  summarise(pGC = median(pGC), depth = median(depth))

gc_min <- min(gc_binned$pGC)
gc_max <- max(gc_binned$pGC)

gc_fit <- mgcv::gam(
  log(depth) ~ s(pGC, k = 5),
  data   =  gc_binned,
  method = "REML"
)

corr_depth <-
  gc_depth %>%
  mutate(
    depth = depth / exp(
      c(predict(
        gc_fit,
        newdata = gc_depth %>% mutate(pGC = pmin(pmax(pGC, gc_min), gc_max)),
        type = "response"  
      ))
    )
  ) %>%
  group_by(chrgrp = chrom == "chrX") %>%
  mutate(
      depth = depth / median(
        depth[depth > 0],
        na.rm = TRUE
      )
  ) %>%
  ungroup()  %>% 
  select(-chrgrp)

WIDTH <- NPHASE * (corr_depth$end[1] - corr_depth$start[1])

phased <-
  corr_depth %>%
  pool_phase_mean(NPHASE = NPHASE) %>%
  filter(end - start == WIDTH)

shards <- parallel::splitIndices(nrow(phased), NSHARDS)

for (i in seq_along(shards)) {
  saveRDS(
    phased$depth[shards[[i]]],
    paste0(output, ".shard_", stringr::str_pad(i, width = nchar(NSHARDS), pad = "0"), ".snorm.rds")
  )

  phased[shards[[i]],] %>%
    select(-depth) %>%
    saveRDS(
      paste0(output, ".shard_", stringr::str_pad(i, width = nchar(NSHARDS), pad = "0"), ".bins.rds")
    )
}



