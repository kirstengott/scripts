#!/usr/bin/Rscript

args <- commandArgs(TRUE)

if (length(args) == 0L || any(c('-h', '--help') %in% args)) {
  message('usage: path/to/busco_plot.R busco_db_path busco_dir
    busco_db_path     path to the busco lineage database 
    busco_dir         path to directory holding busco results to collate in the plot  
    -h, --help        to print help messages')
  q('no')
}

library(tidyverse)
lineage_data <- list.files(args[1], recursive = TRUE, pattern = 'orthogroup_info', full.names = TRUE)

busco_annot <- read_tsv(lineage_data)

data <- list.files(args[2], pattern = 'full', recursive = TRUE, full.names = TRUE) %>%
  purrr::map(., .f = function(x){
    df <- read_delim(x, 
                     comment = "#", 
                     delim = "\t", 
                     col_names = c('Busco', "Status",	"Sequence",	"Score",	"Length"))
    df$name <- basename(x) %>% sub("run_", "", .) %>% 
      sub('full_table_', '', .) %>% 
      sub('.all.maker.proteins.fasta.tsv', '', .) %>% 
      sub('_genomic.faa.tsv', '', .)
    df
  }) %>%
  bind_rows() %>%
  mutate(name = sub(".tsv", "", name))

data_summary <- data %>%
  group_by(name) %>%
  mutate(Total_Busco = length(Status)) %>%
  add_count(Status) %>%
  select(Status, Total_Busco, n) %>%
  distinct() %>%
  mutate(Percent = (n/Total_Busco)*100)


poor_quality <- data_summary %>% filter(Status == 'Complete', Percent < 80) %>% .$name

message(paste('Warning, remove these genomes as they have less than 80% completeness!', 
              paste(poor_quality, collapse = "\n"), sep = '\n'))

ggplot(data_summary, aes(x = name, y = Percent, fill = Status)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = seq(0, 100, by = 10), 
                     expand = c(0,0)) +
  theme_linedraw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggsave(filename = 'busco_genome_completeness.pdf', device = 'pdf', width = 10, height = 20, units = 'in')

data_out <- data_summary %>% select(name, Status, Percent) %>%
  spread(key = Status, value = Percent)

write_csv(data_out, path = "busco_genome_completeness.csv")