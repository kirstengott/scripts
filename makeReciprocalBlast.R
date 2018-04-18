#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/makeReciprocalBlast.R blastdbs outdir blasttype extension num_threads
          blastdbs          path to the blastdb directory
          outdir            path to output directory
          blasttype         type of blast (blastp or blastn)
          extension         fasta extension (fa or fasta)
	  num_threads       cores to use
          -h, --help        to print help messages')
  q('no')
}
library(tools)

args <- c('db', 'blastn' ,'blastn', 'fasta')

blast_type <- args[3]

make_blast <- function(x, y){
  x_mod <- file_path_sans_ext(basename(x))
  y_mod <- file_path_sans_ext(basename(y))
  out   <- paste0(args[2], "/", x_mod, "_", y_mod, ".", blast_type)
  paste(blast_type, '-query', x, '-db', y,
  '-max_target_seqs 1 -evalue 0.001 -num_threads', args[5], '-outfmt 6 -out', out, sep = " ")
}

make_parse <- function(x, y){
  x_mod <- file_path_sans_ext(basename(x))
  y_mod <- file_path_sans_ext(basename(y))
  out   <- paste0(args[2], "/", x_mod, "-", y_mod, ".", blast_type, "Parsed")
  blast_reference <- paste0(args[2], "/", x_mod, "_", y_mod, ".", blast_type)
  blast_invert  <- paste0(args[2], "/", y_mod, "_", x_mod, ".", blast_type)
  paste('parseReciprocalBlast.R', blast_reference, blast_invert, out, sep = " ")
}



fas <- list.files(args[1], full.names = TRUE, pattern = paste0(args[4], '$'))

unique_pairs <- combn(fas, 2, simplify = FALSE)

blast_commands <- unlist(
  lapply(fas, function(x){
  fa_list <- fas[which(!fas == x)]
  make_blast(x, fa_list)
}))

blast_parse_commands <- unlist(
  lapply(unique_pairs, function(x){
    make_parse(x[1], x[2])
}))



write(x = blast_commands, file = 'reciprocal_blasts.sh')
write(x = blast_parse_commands, file = 'parse_reciprocal_blasts.sh')

