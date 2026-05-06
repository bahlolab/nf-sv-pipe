#!/usr/bin/env Rscript

stopifnot(
  require(tidyverse)
)

output   <- commandArgs(trailingOnly = TRUE)[1]
sm_calls <- commandArgs(trailingOnly = TRUE)[-1]
# sm_calls <- list.files('/vast/scratch/users/munro.j/nextflow/work/a0/8dd41ecd4859ae56f9c723b851cb82/', pattern = '.rds', full.names = T)

# merge greedily
MIN_JACC <- 0.9
# recluster if any  pairwise here
JAC_BRK  <- 0.75

samples <- sm_calls  %>% basename() %>% str_remove('\\.(refined_)?calls\\.rds$')

# calculate jaccard index
jaccard <- function(start1, end1, start2, end2) {
  inter <- pmax(0, pmin(end1, end2) - pmax(start1, start2))
  union <- (end1 - start1) + (end2 - start2) - inter
  inter / union
}
# identify sets of calls to merge
merge_calls <- function(calls) {
  calls %>% 
    select(SVTYPE, chrom, start, end, IDX) %>% 
    nest(data = -c(SVTYPE, chrom)) %>% 
    mutate(
      sets = map(data, function(data) {
        jacc <-
          data %>% 
          inner_join(
            data, 
            by = join_by(overlaps(start, end, start, end)),
            suffix = c('_1', '_2')
          ) %>% 
          mutate(jaccard = jaccard(start_1, end_1, start_2, end_2)) %>% 
          select(IDX_1, IDX_2, jaccard)
        
        sets <-
          jacc %>% 
          filter(
            jaccard >= MIN_JACC,
            IDX_1 < IDX_2
          ) %>%
          group_by(IDX_1) %>%
          summarise(set = list(list(unique(sort(c(IDX_1, IDX_2)))))) %>%
          pull(set) %>%
          reduce(function(x, y) {
            ii <- integer()
            for (i in seq_along(x)) {
              if (any(x[[i]] %in% y[[1]])) {
                ii <- c(ii, i)
              }
            }
            if (length(ii)) {
              x[[ii[1]]] <- sort(union(unlist(x[ii]), y[[1]]))
              if (length(ii) > 1) {
                x[ii[-1]] <- NULL
              }
              return(x)
            }
            return(c(x, y))
          }, .init = list(integer())) %>% 
          (function(x) tibble(IDX = x, set = seq_along(x))) %>% 
          unnest(IDX) %>% 
          (function(x) inner_join(
            x, x, by = 'set', 
            suffix = c('_1', '_2'), 
            relationship = 'many-to-many')
          ) %>% 
          inner_join(jacc, by = join_by(IDX_1, IDX_2)) %>% 
          nest(data = -set) %>% 
          mutate(data = map(data, function(data) {
            # split clusters if required
            if (min(data$jaccard) > JAC_BRK) {
              return(
                data %>% 
                  select(IDX = IDX_1) %>% 
                  distinct() %>% 
                  mutate(sub_set = 1L)
              )
            }
            
            hclust <-
              data %>% 
              arrange(IDX_1, IDX_2) %>% 
              mutate(jaccard = 1 - jaccard) %>% 
              pivot_wider(names_from = c(IDX_1), values_from = jaccard) %>% 
              as.data.frame() %>% 
              column_to_rownames('IDX_2') %>% 
              as.dist() %>% 
              hclust(method = 'complete')
            
            for (k in seq(2, length(hclust$labels))) {
              cut <- cutree(hclust, k = k)
              
              jmin <- 
                data %>% 
                mutate(
                  set_1 = cut[as.character(IDX_1)],
                  set_2 = cut[as.character(IDX_2)]
                ) %>% 
                filter(set_1 == set_2) %>% 
                pull(jaccard) %>% 
                min()
              
              if (jmin > JAC_BRK) {
                break
              }
            }
            return(
              data %>% 
                select(IDX = IDX_1) %>% 
                distinct() %>% 
                mutate(sub_set = cut[as.character(IDX)])
            )
          })
          ) %>% 
          unnest(data) %>%
          select(-any_of('data')) %>% 
          bind_rows(tibble(sub_set = integer())) %>%
          mutate(sub_set = replace_na(sub_set, 1L)) %>% 
          mutate(set = str_c(set, sub_set, sep = '-'))
      })
    ) %>% 
    unnest(sets) %>% 
    mutate(set = str_c(chrom, SVTYPE, set, sep = '-')) %>% 
    select(IDX, set)
}

CALLS <-
  map_df(sm_calls, readRDS) %>%
  # mutate(precise = precise_left & precise_right) %>% 
  mutate(precise = FALSE) %>% 
  mutate(
    SVLEN = end - start,
    IDX = row_number()
  )

MERGE_SETS_1 <-
  CALLS %>% 
  merge_calls()

# exclude excessive calls from outlier samples
ncall_threshold <-
  CALLS %>% 
  count(sample, SVTYPE) %>% 
  full_join(expand_grid(SVTYPE = unique(CALLS$SVTYPE), sample=samples)) %>% 
  mutate(
    n = replace_na(n, 0),
    x = log1p(n)
  ) %>% 
  group_by(SVTYPE) %>% 
  summarise(
    threshold = median(x, na.rm = TRUE) + -qnorm(0.01) * mad(x, na.rm = TRUE),
    threshold = ceiling(exp(threshold) - 1),
    .groups = 'drop'
  )

# remove excess calls from outlier samples
# use model based on length/precision/CN/nmerged
FLT_CALLS <-
  CALLS %>% 
  mutate(SVLEN = end - start) %>% 
  select(sample, SVTYPE, SVLEN, CN, SNR, precise, IDX) %>% 
  left_join(
    MERGE_SETS_1 %>%
      add_count(set, name = 'nmerged'),
    by = 'IDX'
  ) %>% 
  mutate(nmerged = replace_na(nmerged, 1)) %>% 
  inner_join(ncall_threshold) %>% 
  add_count(sample, SVTYPE) %>% 
  mutate(is_outlier = n > threshold) %>% 
  nest(data = -SVTYPE) %>% 
  mutate(
    data = map(data, function(data) {
      if (sum(data$is_outlier) == 0) {
        return(mutate(data, p_outlier = 0))
      }
      fit <- glm(
        is_outlier ~ log(SVLEN) + precise + CN + SNR + log(nmerged),
        data = data,
        family = binomial()
      )
      mutate(
        data,
        p_outlier = predict(fit, type = 'response')
      )
    })
  ) %>% 
  unnest(data) %>% 
  arrange(sample, SVTYPE, p_outlier) %>% 
  group_by(sample, SVTYPE) %>% 
  filter(row_number() < threshold) %>% 
  ungroup()

MERGE_SETS_2 <-
  CALLS %>% 
  semi_join(FLT_CALLS, by = 'IDX') %>% 
  merge_calls()

# find calls to merge
MERGED_CALLS_LONG <-
  bind_rows(
    CALLS %>% 
      semi_join(FLT_CALLS, by = 'IDX') %>% 
      anti_join(MERGE_SETS_2, by = 'IDX'),
    CALLS %>% 
      inner_join(MERGE_SETS_2, by = 'IDX') %>% 
      group_by(set) %>% 
      mutate(
        start = round(weighted.mean(start, if_else(precise, 2, 1))),
        end   = round(weighted.mean(end,   if_else(precise, 2, 1))),
        precise = sum(precise) / n() > 0.80
      ) %>% 
      ungroup()
  ) %>% 
  mutate(SVLEN = end - start) %>% 
  select(chrom, start, end, SVTYPE, SVLEN, precise, sample, GT, CN) %>% 
  arrange(chrom, start, end, SVTYPE) %>% 
  group_by(chrom, start, end, SVTYPE) %>% 
  mutate(AC = sum(str_count(GT, '1')),
         AN = 2 * length(samples),
         AF = AC / AN,
         .after = precise
  ) %>% 
  ungroup() %>% 
  group_by(chrom, start, end, SVTYPE, SVLEN, precise, AC, AN, AF, sample) %>% 
  slice(1) %>% 
  ungroup()

saveRDS(MERGED_CALLS_LONG, 'merged_calls.rds')

# write VCF
VCF <-
  MERGED_CALLS_LONG %>% 
  bind_rows(tibble(sample = samples, dummy = TRUE)) %>% 
  mutate(format = str_c(GT, ':', round(CN, 2))) %>% 
  select(-GT, -CN) %>% 
  pivot_wider(
    names_from = sample,
    values_from = format,
    values_fill = '0/0:.',
  ) %>% 
  filter(is.na(dummy)) %>% 
  (function(x) {
    transmute(
      x,
      CHROM  = chrom,
      POS    = as.integer(start) + 1L,
      ID     = '.',
      REF    = "N",
      ALT    = str_c("<", SVTYPE, ">"),
      QUAL   = ".",
      FILTER = "PASS",
      INFO   = str_c(
        "SVTYPE=", SVTYPE, 
        ";END=", end,
        ";SVLEN=", if_else(SVTYPE == "DEL", -SVLEN, SVLEN),
        ';PRECISE=', as.integer(precise)
      ),
      FORMAT = "GT:CN",
    ) %>% 
      bind_cols(select(x, all_of(samples)))
  })

# header
c(
  "##fileformat=VCFv4.2",
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
  '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
  '##INFO=<ID=PRECISE,Number=1,Type=Integer,Description="Breakpoints are precise">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number estimate">',
  str_c("#", str_c(colnames(VCF), collapse = "\t"))
) %>% 
  write_lines(output)

write_tsv(VCF, output, col_names = F, append = TRUE)

