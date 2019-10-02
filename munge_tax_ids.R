#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/munge_tax_ids.R metadata_tax seqid_map
          metadata_tax      taxonomy metadata file
          seqid_map         map of unique tree ids to sequence ids
          -h, --help        to print help messages')
  q('no')
}


library(tidyverse)


munge_accessions <- function(x){
  paste(strsplit(x, split = "_")[[1]][c(1,2)], collapse = "_")
}

metadata <- read_tsv(args[1], col_names = FALSE) %>%
  select(X1, X28, X29, X30) %>%
  rename('acc' = X1,
         'family' = X28,
         'genus' = X29,
         'genus_species' = X30) %>%
  distinct() %>%
  mutate(genus_species = sub("\\]", "", sub("\\[", "", genus_species)))


ids <- read_tsv(args[2], col_names = c('genome_id', 'id')) %>%
  group_by(genome_id) %>%
  mutate(acc = munge_accessions(genome_id))


ids %>% left_join(., metadata, by = 'acc') %>%
  distinct() %>%
  write_tsv(., path = 'id_metadata.txt')

