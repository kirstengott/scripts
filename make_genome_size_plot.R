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


find_seq_error <- function(vector_nmer){
  ## this function takes a vector of nmers and finds the point at which it reverses
  ## and returns the reversal nmer value
  n <- vector_nmer[1]
  for(y in vector_nmer){
    if(y <= n){
      n = y
    } else {
      reverse = y
      return(reverse)
      break
    }
  }
}

# find_seq_error2 <- function(vector1, vector2){
#   n <- vector1[1]
#   for(y in vector1){
#     if(y <= n){
#       n = y
#     } else {
#       reverse = y
#       return(vector2[which(vector1 == reverse)])
#       break
#     }
#   }
# }

library(tidyverse, verbose = FALSE, quietly = TRUE)

#args <- c('jellyfish/new_output/PE1_31mer.histo', 'test', 'pe3')
#args <- c("jellyfish/hist/1054.histo", 'test', '1054')

x <- read_delim(args[1], col_names = c('Coverage', 'Nmers'), delim = ' ') %>%
                    mutate(name = sub(".histo", "", basename(args[1])))


plot_title <- basename(args[1])


## find where the sequenceing error ends
seq_error_end <- find_seq_error(x$Nmers)

index_error_end <- which(x$Nmers == seq_error_end)

## return the index of where the sequencing error ends
coverage_start <- x[x$Nmers == seq_error_end, 'Coverage', drop = TRUE]

## filter the datatable for positions that do not have sequencing error
x_filt <- x %>%
  filter(Coverage >= coverage_start,
         Coverage <= nrow(x))

## find the maximum height of the single copy region
peak_cov <- x_filt %>% 
  filter(Nmers == max(Nmers)) %>% 
  .$Coverage



## pull out the maximum value of the y axis and arbitrarily add 1000 for plotting purposes
y_max    <- (x_filt %>% filter(Nmers == max(Nmers)) %>% .$Nmers) + 1000


## size of the single copy region of the genome
single_copy_region <- round((sum(x_filt$Coverage*x_filt$Nmers)/peak_cov)/1000000, digits = 2)

## total genome size
total_genome_size <- round((sum(as.numeric(x[index_error_end:nrow(x),1, drop = TRUE]*x[index_error_end:nrow(x),2, drop = TRUE]))/peak_cov)/1000000, digits = 2)



ggplot(x, aes(x = Coverage, y = Nmers)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 200, by = 10), labels = seq(0, 200, by = 10),
                     limits = c(0, 200))+
  ylim(0, y_max) +
  theme_kirsten(angle=0) +
  labs(title = plot_title,
       subtitle = paste('Estimated Genome Size:', total_genome_size),
       caption = paste('Peak Coverage:', peak_cov)) +
  ggsave(device = 'pdf', filename = paste0(basename(args[1]), ".pdf"))




write(x = paste(basename(args[1]), total_genome_size, single_copy_region, sep = "\t"), file = 'genome_sizes.txt', append = TRUE)
