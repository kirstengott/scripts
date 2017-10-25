#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/makeReciprocalBlast.R blastdbs outdir
          blastdbs          path to the blastdbs
          outdir            path to output
          blasttype         type of blast
          -h, --help        to print help messages')
  q('no')
}
library(tools)

#args <- c('../data/300_genomes_proteins/', 'recips', 'blastp')  

blast_type <- args[3]

make_blast <- function(x, y){
  x_mod <- file_path_sans_ext(basename(x))
  y_mod <- file_path_sans_ext(basename(y))
  out   <- paste0(args[2], "/", x_mod, "_", y_mod, ".", blast_type)
  paste(blast_type, '-query', x, '-db', y, 
  '-max_target_seqs 1 -evalue 0.001 -num_threads 6 -outfmt 6 -out', out, sep = " ")
}

make_parse <- function(x, y){
  x_mod <- file_path_sans_ext(basename(x))
  y_mod <- file_path_sans_ext(basename(y))
  out   <- paste0(args[2], "/", x_mod, "_", y_mod, ".", blast_type, "_parsed")
  blast_reference <- paste0(args[2], "/", x_mod, "_", y_mod, ".", blast_type)
  blast_invert  <- paste0(args[2], "/", y_mod, "_", x_mod, ".", blast_type)
  paste('parseReciprocalBlast.R', blast_reference, blast_invert, sep = " ")
}

fas <- list.files(args[1], full.names = TRUE, pattern = 'fas$')

blast_commands <- unlist(
  lapply(fas, function(x){
  fa_list <- fas[which(!fas == x)]
  make_blast(x, fa_list)
}))

parse_commands <- unlist(
  lapply(fas, function(x){
    fa_list <- fas[which(!fas == x)]
    make_parse(x, fa_list)
}))




blast_parse_commands <- unlist(
  lapply(fas, function(x){
    fa_list <- fas[which(!fas == x)]
    make_blast(x, fa_list)
  }))


write(x = blast_commands, file = 'reciprocal_blasts.txt')
write(x = blast_parse_commands, file = 'parse_reciprocal_blasts.txt')

