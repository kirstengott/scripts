#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
    message('usage: path/to/parseReciprocalBlast.R referenceBlast invertBlast output
    referenceBlast             Reference blast output.
    invertBlast                Reciprocal of reference blast output.
    output                     Name of outfile.
    -h, --help        to print help messages')
    q('no')
}

library(dplyr)
library(tidyr)


reference <- data.table::fread(args[1], data.table = FALSE)
invert    <- data.table::fread(args[2], data.table = FALSE)


cols <- c("query", 
          "subject", 
          "perc_id", 
          "align_len", 
          "mismatches", 
          "gaps", 
          "q_start", 
          "q_end", 
          "s_start", 
          "s_end", 
          "e_val", 
          "bit_score")

colnames(reference) <- paste0('reference_', cols)
colnames(invert)    <- paste0('invert_', cols)



referenceRank <- reference %>% 
  filter(reference_e_val <= 1e-3) %>% 
  group_by(reference_query) %>% 
  mutate(rank1 = floor(rank(-reference_bit_score))) %>%
  ungroup() %>% 
  filter(rank1 == 1) %>% 
  rename('common_subject'= reference_subject)

invertRank <- invert %>% 
  filter(invert_e_val <= 1e-3) %>% 
  group_by(invert_query) %>% 
  mutate(rank2 = floor(rank(-invert_bit_score))) %>%
  ungroup()%>% 
  rename('common_subject' = invert_query, 
         'reference_query' = invert_subject)


ranker    <- left_join(referenceRank, invertRank)
rankfinal <- ranker %>% select(reference_query, 
                               common_subject, 
                               reference_e_val, 
                               rank2)



write.table(na.omit(rankfinal), 
            file = args[3], 
            quote = FALSE, 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE)
