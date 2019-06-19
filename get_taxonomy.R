#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
    message('usage: path/to/get_taxonomy input output
    input             tsv with a taxid column to indicate ncbi taxons
    output            Name of outfile
    -h, --help        to print help messages')
    q('no')
}


library('taxonomizr')
library(tidyverse)

metadata <- read_tsv(args[1])

taxid <- metadata %>% .$taxid %>% unique()

getNamesAndNodes()
getAccession2taxid()
read.accession2taxid(list.files('.','accession2taxid.gz$'),'accessionTaxa.sql')
taxaNodes <-read.nodes('nodes.dmp')
taxaNames <-read.names('names.dmp')

taxa_df <- getTaxonomy(taxid,taxaNodes,taxaNames)

taxa_df <- taxa_df %>% data.frame(., stringsAsFactors = FALSE) %>%
  rownames_to_column(var = 'taxid') %>%
  mutate(taxid = as.integer(taxid))


all_df <- left_join(metadata, taxa_df, by = 'taxid')

write_tsv(all_df, path = args[2])


