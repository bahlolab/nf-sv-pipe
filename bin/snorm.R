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

depth <-
  read_tsv(
    input,
    col_names = c("chrom", "start", "end", "depth"),
    col_types = "ciin"
  ) %>%
  filter(chrom %in% paste0("chr", c(1:22, "X"))) %>%
  group_by(chrom == "chrX") %>%
  mutate(
    depth = depth / median(
      depth[depth > 0],
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  pool_phase_mean(NPHASE = NPHASE)

shards <- parallel::splitIndices(nrow(depth), NSHARDS)

for (i in seq_along(shards)) {
  saveRDS(
    depth$depth[shards[[i]]],
    paste0(output, ".shard_", stringr::str_pad(i, width = nchar(NSHARDS), pad = "0"), ".snorm.rds")
  )

  depth[shards[[i]],] %>%
    select(-depth) %>%
    saveRDS(
      paste0(output, ".shard_", stringr::str_pad(i, width = nchar(NSHARDS), pad = "0"), ".bins.rds")
    )
}



