#!/usr/bin/env Rscript

stopifnot(
    require(readr),
    require(dplyr),
    require(tidyr),
    require(stringr),
    require(purrr)
)

output     <- commandArgs(trailingOnly = TRUE)[1]
calls      <- commandArgs(trailingOnly = TRUE)[2]
bpt_depth  <- commandArgs(trailingOnly = TRUE)[3]
WINDOW     <- as.integer(commandArgs(trailingOnly = TRUE)[4])
WIDTH      <- as.integer(commandArgs(trailingOnly = TRUE)[5])

WINDOW <- 2*floor(WINDOW/2)
MAX_DEL <- 1.5
MIN_DUP <- 2.5

CALLS <- 
  readRDS(calls) %>% 
  mutate(IDX = row_number())

if (nrow(CALLS) == 0) {
    saveRDS(CALLS, paste0(output, ".refined_calls.rds"))
    q("no", status = 0)
}

CALL_DP <-
  CALLS %>% 
  rowwise() %>% 
  reframe(
    bind_rows(
      tibble(
        IDX,
        SVTYPE, 
        chrom, 
        BKPT = 'LEFT',
        pos = seq(start - WINDOW/2, length.out = WINDOW),
        rel_pos = pos - start
      ),
      tibble(
        IDX,
        SVTYPE, 
        chrom, 
        BKPT = 'RIGHT',
        pos = seq(end - WINDOW/2, length.out = WINDOW),
        rel_pos = pos - end
      )
    )
  ) %>% 
  left_join(
    read_tsv(
      bpt_depth,
      col_names = c("chrom", "pos", "depth"),
      col_types = "cii"
    ) %>% 
      distinct(), 
  by = join_by(chrom, pos))

CALL_SMRY <-
  CALL_DP %>% 
  group_by(IDX, SVTYPE, BKPT) %>% 
  summarise(
    DPOUT =  median(depth[(rel_pos < 0 & BKPT == "LEFT") | (rel_pos > 0 & BKPT == "RIGHT")], na.rm = T),
    DPIN  = median(depth[!((rel_pos < 0 & BKPT == "LEFT") | (rel_pos > 0 & BKPT == "RIGHT"))], na.rm = T),
    .groups = 'drop'
  ) %>% 
  filter(
    (SVTYPE == 'DEL' & 2*DPIN/DPOUT < MAX_DEL) |
      (SVTYPE == 'DUP' & 2*DPIN/DPOUT > MIN_DUP),
    DPOUT > 0
  )

CALL_DP_FLT <-
  CALL_DP %>% 
  inner_join(CALL_SMRY) %>%  
  mutate(IDX2 = row_number()) 


FINE <-
  CALL_DP_FLT %>% 
  mutate(depth = replace_na(depth, 0L)) %>% 
  mutate(
    med_left   = data.table::frollmedian(depth, WIDTH, na.rm = T),
    mean_left  = data.table::frollmean(depth, WIDTH, na.rm = T),
    med_right  = lead(data.table::frollmedian(depth, WIDTH, na.rm = T), WIDTH),
    mean_right = lead(data.table::frollmean(depth, WIDTH, na.rm = T)  , WIDTH),
  ) %>% 
  filter(
    rel_pos >= -WINDOW/2 + WIDTH,
    rel_pos <   WINDOW/2 - WIDTH + 1
  ) %>% 
  mutate(
    RELMED = if_else(
      BKPT == 'LEFT',
      2 * med_right / med_left,
      2 * med_left / med_right,
    ),
    RELMEAN = if_else(
      BKPT == 'LEFT',
      2 * mean_right / mean_left,
      2 * mean_left / mean_right,
    )
  ) %>% 
  filter(
    (SVTYPE == 'DEL' & RELMED < MAX_DEL & RELMEAN < MAX_DEL) |
      (SVTYPE == 'DUP' & RELMED > MIN_DUP & RELMEAN > MIN_DUP),
    DPOUT > 0
  )  %>% 
  group_by(IDX, SVTYPE, BKPT) %>% 
  filter(
    (SVTYPE == 'DEL' & RELMED == min(RELMED)) |
      (SVTYPE == 'DUP' & RELMED == max(RELMED))
  ) %>% 
  filter(
    (SVTYPE == 'DEL' & RELMEAN == min(RELMEAN)) |
      (SVTYPE == 'DUP' & RELMEAN == max(RELMEAN))
  ) %>% 
  group_by(IDX, IDX2, SVTYPE, BKPT, pos) %>% 
  summarise(
    res = {
      # rank sum statistic
      rank <- CALL_DP_FLT$depth[seq(IDX2-WIDTH, length.out = WIDTH*2)]
      abs(
        sum(rank[seq(1, length.out = WIDTH)]) - 
          sum(rank[seq(1+WIDTH, length.out = WIDTH)])
      )
    },
    .groups = 'drop'
  ) %>% 
  group_by(IDX, BKPT) %>%
  slice_min(order_by = res, with_ties = F) %>% 
  ungroup()

CALLS %>% 
  left_join(
    FINE %>% 
      ungroup() %>% 
      filter(BKPT == 'LEFT') %>% 
      transmute(IDX, new_start = pos)
  ) %>% 
  left_join(
    FINE %>% 
      ungroup() %>% 
      filter(BKPT == 'RIGHT') %>% 
      transmute(IDX, new_end = pos)
  ) %>% 
  mutate(
    precise_left = !is.na(new_start),
    precise_right= !is.na(new_end),
    start = if_else(!is.na(new_start), new_start, start),
    end   = if_else(!is.na(new_end), new_end, end),
    SVLEN = end - start,
  ) %>% 
  select(
    -IDX, -new_start, -new_end
  ) %>% 
  saveRDS(paste0(output, ".refined_calls.rds"))

