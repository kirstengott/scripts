#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/split_ragoo.R input output
          input             Ragoo output.
          output            Name of outfile
          -h, --help        to print help messages')
  q('no')
}

library(Biostrings)
library(tidyverse)



ragoo <- readDNAStringSet(args[1])

t <- Biostrings::strsplit(x = ragoo, split = "N")

chr1_split <- t$Chr0_RaGOO
chr_rest   <- ragoo[-1]


chr_all <- c(chr_rest, chr1_split)
names(chr_all) <- seq(1, length(chr_all))


sort_df <- data.frame(names = names(chr_all), width = width(chr_all)) %>% arrange(-width)

chr_sort        <- chr_all[sort_df$names]
names(chr_sort) <- paste0('node_', seq(1, length(chr_sort)), "_length_", width(chr_sort))


### remove sequences with only Ns
all_n_seq <- sapply(seq(1, length(chr_sort)), function(y){
  x <- chr_sort[y]
  t <- alphabetFrequency(x)
  if(t[,'N'] == width(x)){
    return(y)
  } else {
    return(NA)
    }
})


all_ns <- na.omit(all_n_seq)

if(length(all_ns) >0){
  chr_final <- chr_sort[-all_ns]
} else {
  chr_final <- chr_sort
}


writeXStringSet(chr_final, args[2])


