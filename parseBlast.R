#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
    message('usage: path/to/parseBlast.R input output
    input             Blast output.
    output            Name of outfile
    -h, --help        to print help messages')
    q('no')
}

library(dplyr)


blast           <- data.table::fread(args[1], sep="\t", stringsAsFactors=FALSE, data.table = FALSE)

colnames(blast) <- c("query", "subject", "perc_id", "align_len", "mismatches", "gaps", "q_start", "q_end", "s_start", "s_end", "e_val", "bit_score")

options(dplyr.width = Inf)


blast.final <- blast %>% filter(e_val< 0.001) %>% 
    group_by(query) %>% 
        filter(rank(-bit_score, ties.method="first") == 1) %>%
            ungroup()


if (nrow(blast.final) == 0) {
    message("No evalues < 0.001.")
} else {
## write out the results
write.table(blast.final[ ,c("query", "subject", "e_val", "perc_id", "align_len")], file = args[2], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

