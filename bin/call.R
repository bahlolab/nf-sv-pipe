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

MIN_SEED <- 3L # CNVs grown out from seeds of a least 3 adjacent bins
MIN_CALL <- 5L # number of bins required to keep a call
MAX_GAP <- 0.10
MIN_P   <- 0.001
MAX_DEL <- 1.5
MIN_DUP <- 2.5
MAX_DUP <- 5

# hard limit on number of class - most samples should be unaffected
MAX_N_DEL <- 250
MAX_N_DUP <- 250

bnorm <- map_df(inputs, readRDS) %>% arrange(phase, chrom, start)
width <- bnorm$end[1] - bnorm$start[1]

status <-
  bnorm %>%
  filter(
    is.finite(Z),
  ) %>%
  mutate(
    IDX = row_number(),
    B_PVAL = 2 * pnorm(-abs(Z)),
    ZZ = (Z - median(Z, na.rm = TRUE)) / mad(Z, na.rm = TRUE),
    S_PVAL = 2 * pnorm(-abs(ZZ))
  ) %>%
  ungroup() %>%
  mutate(SVTYPE = case_when(
    B_PVAL > MIN_P ~ NA_character_,
    S_PVAL > MIN_P ~ NA_character_,
    CN < MAX_DEL   ~ "DEL",
    (CN > MIN_DUP & CN < MAX_DUP) ~ "DUP",
    TRUE           ~ NA_character_
  ))

seeds <-
  status %>%
  filter(!is.na(SVTYPE)) %>%
  mutate(
    BREAK = SVTYPE != lag(SVTYPE) | start != lag(end) | round(CN) != lag(round(CN)),
    SEG = cumsum(replace_na(BREAK, TRUE)),
  ) %>%
  group_by(SEG) %>%
  filter(n() >= MIN_SEED) %>%
  summarise(
    chrom = first(chrom),
    phase = first(phase),
    SVTYPE = first(SVTYPE),
    start = first(start),
    end = last(end),
    med_CN = median(CN),
    CN = list(CN),
    n = n(),
    .groups = 'drop'
  )

# merge adjacent SEED regions to give calls per phase
phase_calls <-
  seeds %>%
  nest(data = -c(chrom, phase)) %>%
  mutate(data = map(data, function(data) {
    stage_2 <- FALSE
    while (TRUE) {
      n <- nrow(data)
      if (n == 0) { break }
      data <-
        data %>%
        mutate(
          GAP = (start - lag(end)) / (end - lag(start)),
          BREAK = SVTYPE != lag(SVTYPE) | abs(med_CN - lag(med_CN)) > 0.25 | GAP > MAX_GAP,
          SEG = cumsum(replace_na(BREAK, TRUE)),
        ) %>%
        group_by(SEG) %>%
        summarise(
          SVTYPE = first(SVTYPE),
          start = first(start),
          end = last(end),
          CN = list(unlist(CN)),
          med_CN = median(unlist(CN)),
          r_CN = round(med_CN),
          n = lengths(CN),
          .groups = "drop"
        )
      if (nrow(data) == n) {
        # break
        if (stage_2) { break }
        stage_2 <- TRUE
        data <- filter(data, n >= MIN_CALL)
      }
    }
    return(data)
  })) %>%
  unnest(data) %>% 
  select(-SEG)

# merge calls across phases greedily
merged <-
  phase_calls %>% 
  select(chrom, SVTYPE, r_CN, start, end) %>% 
  (function(x) {
    while(TRUE) {
      n <- nrow(x)
      x <- 
        mutate(x, index = row_number()) %>% 
        left_join(
          x %>% 
            select(chrom, SVTYPE, start, end, r_CN),
          by = join_by(chrom, SVTYPE, r_CN, overlaps(start, end, start, end)),
          suffix = c('', '_2')
        ) %>% 
        group_by(chrom, SVTYPE, r_CN, index) %>% 
        summarise(
          start = min(start_2),
          end = max(end_2),
          .groups = 'drop'
        ) %>% 
        select(-index) %>% 
        distinct()
      if (nrow(x) == n) { return (x) }
    }
  }) %>% 
  left_join(
    phase_calls,
    by = join_by(chrom, SVTYPE, r_CN, overlaps(start, end, start, end)),
    suffix = c('', '_2')
  ) %>% 
  group_by(chrom, SVTYPE, start, end, r_CN) %>% 
  summarise(
    med_CN = median(unlist(CN)),
    .groups = 'drop'
  ) %>% 
  mutate(IDX = row_number())

# Trim outlier bins from ends for better accuracy
trimmed <- 
  merged %>% 
  inner_join(
    status %>% 
      select(chrom, start, end, CN),
    by = join_by(chrom, overlaps(start, end, start, end)),
    suffix = c('', '_2')
  ) %>% 
  filter(start_2 >= start, end_2 <= end) %>% 
  arrange(chrom, start_2, end_2) %>% 
  group_by(chrom, SVTYPE, start, end) %>% 
  (function(data) {
    while(TRUE) {
      n <- nrow(data)
      data <- 
        data %>% 
        mutate(
          index = row_number(),
          Z = (CN - median(CN, na.rm = TRUE)) / mad(CN, na.rm = TRUE)
        ) %>% 
        filter(
           ! ((index == 1 | index == n()) &
                abs(Z) > -qnorm(0.01) & 
                is.finite(Z) &
                (SVTYPE[1] == 'DEL' & Z > 0 | SVTYPE[1] == 'DUP' & Z < 0)
           )
        ) %>% 
        mutate(start = min(start_2), end = max(end_2))
      if (nrow(data) == n) { return (data) }
    }
  }) %>% 
  summarise(
    med_CN = median(CN),
    r_CN = round(med_CN),
    .groups = 'drop'
  )

snr <-
  trimmed %>% 
  mutate(svlen = end - start,
         s = start - svlen %/% 2,
         e  = end + svlen %/% 2
  ) %>% 
  inner_join(
    status %>% 
      select(chrom, start, end, CN),
    by = join_by(chrom, overlaps(s, e, start, end)),
    suffix = c('', '_2')
  ) %>% 
  mutate(
    status = case_when(
      end_2 <= start                                  ~ 'left_out',
      start_2 >= start & end_2 < start + svlen %/% 2  ~ 'left_in',
      start_2 >= start + svlen %/% 2 & end_2 <= end ~   'right_in',
      start_2 >= end                                  ~ 'right_out'
    )
  ) %>% 
  filter(!is.na(status)) %>% 
  group_by(
    chrom, SVTYPE, start, end, med_CN, r_CN, status
  ) %>% 
  summarise(
    mu = mean(CN),
    sigma = sd(CN),
    n = n(),
    .groups = 'drop'
  ) %>% 
  pivot_wider(names_from = status, values_from = c(mu, sigma, n)) %>% 
  bind_rows(tibble(
    mu_right_in = numeric(),
    mu_right_out = numeric(),
    mu_left_in = numeric(),
    mu_left_out = numeric(),
    sigma_right_in = numeric(),
    sigma_right_out = numeric(),
    sigma_left_in = numeric(),
    sigma_left_out = numeric(),
    n_right_in = numeric(),
    n_right_out = numeric(),
    n_left_in = numeric(),
    n_left_out = numeric(),
    )
  ) %>% 
  mutate(
    SNR_right = abs(mu_right_in - mu_right_out) / 
      sqrt(sigma_right_in**2 / n_right_in + sigma_right_out**2 / n_right_out),
    SNR_left = abs(mu_left_in - mu_left_out) / 
      sqrt(sigma_left_in**2 / n_left_in + sigma_left_out**2 / n_left_out),
    SNR = (coalesce(SNR_right, SNR_left) + coalesce(SNR_left, SNR_right)) / 2
    ) %>% 
  select(chrom, SVTYPE, start, end, med_CN, r_CN, SNR, SNR_right, SNR_left) %>% 
  filter(!is.na(SNR))

calls <-
  snr %>% 
  group_by(SVTYPE) %>% 
  arrange(desc(SNR)) %>% 
  filter(row_number() <= if_else(SVTYPE == 'DEL', MAX_N_DEL, MAX_N_DUP)) %>% 
  ungroup() %>% 
  transmute(
    sample = output,
    chrom,
    start,
    end,
    SVTYPE,
    SVLEN = end - start,
    CN = med_CN,
    SNR = SNR,
    GT  = case_when(
      r_CN == 0 ~ '1/1',
      r_CN == 1 ~ '0/1',
      r_CN == 3 ~ '0/1',
      r_CN >= 4 ~ '1/1',
      TRUE      ~ '0/0'
    )
  ) %>% 
  arrange(chrom, start, end, SVTYPE)


saveRDS(calls, paste0(output, ".calls.rds"))

breakpoint_windows <-
  bind_rows(
    calls %>%
      transmute(chrom, end = start + width %/% 2, start = start - width %/% 2),
    calls %>%
      transmute(chrom, start = end - width %/% 2, end = end + width %/% 2),
  ) %>%
  select(chrom, start, end) %>%
  arrange_all() %>%
  transmute(reg = paste0(chrom, ":", start, "-", end))

write_tsv(breakpoint_windows, paste0(output, ".bpt.txt"), col_names = FALSE)


