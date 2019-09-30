#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/make_genome_size_plot.R jellyfish.hist output_prefix title
          jellyfish.hist    output kmer histogram from jellyfish
          output_prefix     outfile prefix
          title             plot title
          -h, --help        to print help messages')
  q('no')
}

source('~/scripts/theme_kirsten.R')

find_seq_error <- function(vector){
  n <- vector[1]
  for(y in vector){
    if(y <= n){
      n = y
    } else {
      reverse = y
      return(reverse)
      break
    }
  }
}


find_seq_error2 <- function(vector1, vector2){
  n <- vector1[1]
  for(y in vector1){
    if(y <= n){
      n = y
    } else {
      reverse = y
      return(vector2[which(vector1 == reverse)])
      break
    }
  }
}

library(tidyverse, verbose = FALSE, quietly = TRUE)

#args <- c('jellyfish/new_output/PE1_31mer.histo', 'test', 'pe3')


args <- c("1054_histo_orig")

x <- read_delim(args[1], col_names = c('Coverage', 'Nmers'), delim = ' ') %>%
                    mutate(name = sub(".histo", "", basename(args[1])))


plot_title <- args[3]



seq_error_end <- find_seq_error(x$Nmers)
coverage_start <- x[x$Nmers == seq_error_end, 'Coverage', drop = TRUE]


x_filt <- x %>%
  filter(Coverage >= coverage_start,
         Coverage <= nrow(x))


peak_cov <- x_filt %>% filter(Nmers == max(Nmers)) %>% .$Coverage
y_max    <- (x_filt %>% filter(Nmers == max(Nmers)) %>% .$Nmers) + 1000


## number of k-mers in the histogram
genome_size <- round((sum(x_filt$Coverage*x_filt$Nmers)/peak_cov)/1000000, digits = 2)


ggplot(x, aes(x = Coverage, y = Nmers)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 200, by = 10), labels = seq(0, 200, by = 10),
                     limits = c(0, 200))+
  ylim(0, y_max) +
  theme_kirsten(angle=0) +
  labs(title = plot_title,
       subtitle = paste('Estimated Genome Size:', genome_size),
       caption = paste('Peak Coverage:', peak_cov)) +
  ggsave(device = 'pdf', filename = paste0(args[2], ".pdf"))




write(x = paste(basename(args[1]), genome_size, sep = "\t"), file = 'genome_sizes.txt', append = TRUE)
