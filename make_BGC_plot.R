library(ggplot2)
library(gggenes)
library(GenomicRanges)


## colors pulles from the stata pallete 's2' in ggthemes package
color_mapping <- c('ACP' =  "#1a476f",
                   'AMP-binding' = "#90353b",
                   'Aminotran_1_2' = "#55752f",
                   'Condensation'= "#e37e00",
                   'ECH' = "#6e8e84",
                   'MT' = "#c10534",
                   'NAD_binding_4' = "#938dd2",
                   'PCP' = "#cac27e",
                   'AT' = "#a0522d",
                   'DH' = "#7b92a8",
                   'ER' = "#2d6d66",
                   'KR' = "#9c8847",
                   'KS' = "#bfa19c",
                   'TD' = "#ffd200",
                   'Thioesterase' = "#d9e6eb",
                   'Additional Biosynthetic Domains' = "#404040")

# color_mapping2 <- c('ACP' =  "#1a476f",
#                     'AMP-binding' = "#90353b",
#                     'Aminotran_1_2' = "#55752f",
#                     'Condensation'= "#e37e00",
#                     'ECH' = "#6e8e84",
#                     'MT' = "#c10534",
#                     'NAD_binding_4' = "#938dd2",
#                     'PCP' = "#cac27e",
#                     'PKS_AT' = "#a0522d",
#                     'PKS_DH' = "#7b92a8",
#                     'PKS_ER' = "#2d6d66",
#                     'PKS_KR' = "#9c8847",
#                     'PKS_KS' = "#bfa19c",
#                     'TD' = "#ffd200",
#                     'Thioesterase' = "#d9e6eb",
#                     'Additional Biosynthetic Domains' = "#404040")
# 
#               

#gff3_file = outfile
make_BGC_plot <- function(gff3_file, genome_levels){

  genome_levels <- meta_levels
  
  # lapply(gff$Note, function(x){
  #   sub("GenomeID:", "", rev(x)[1])
  # })
  
  
  
  gff <- rtracklayer::import(gff3_file, format = 'gff3') %>%
    data.frame() %>%
    mutate(direction = ifelse(strand == "+", yes = 1, no = 
                                ifelse(strand == "-", yes = -1, no = NA))) %>%
    mutate(seqnames = sub(" .*$", "", seqnames)) %>%
    group_by(seqnames, start, end) %>%
    mutate(genome = sub("GenomeID:", "", Note[[1]][which(grepl('GenomeID:', Note[[1]]))])) %>%
    ungroup()
  
  gff$genome <- factor(gff$genome, levels = genome_levels)
  
  product <- sub("product:", "", gff$Note[[1]][which(grepl('product:', gff$Note[[1]]))])

  ## may integrate later somehow, not sure how
  core_biosynthetic <- gff %>% filter(type == 'aSDomain') %>% 
    mutate(ID = sub("PKS_", "", ID)) %>% GRanges()
  
  ## may need to make more consice later
  subdomains <- gff %>% filter(type == 'PFAM_domain') %>% GRanges()
  
  overlaps <- mergeByOverlaps(query = subdomains, subject = core_biosynthetic) %>% data.frame()
  
  
sd_data <-  gff %>% filter(type == 'PFAM_domain') %>% 
    mutate(color = ifelse(ID %in% overlaps$subdomains.ID, 
           yes = overlaps$core_biosynthetic.ID, 
           no = 'Additional Biosynthetic Domains')) %>%
  mutate(color = ifelse(ID %in% unique(core_biosynthetic$ID),
                        yes = ID,
                        no = color))
  
  #print(unique(sd_data$color))
  gff %>% 
    filter(type %in% c('CDS')) %>%
    ggplot(aes(xmin = start, 
               xmax = end, 
               y = genome, 
               forward = direction)) +
    geom_gene_arrow(fill = 'white',
                    arrowhead_height = unit(3, "mm"), 
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(3, 'mm'),
                    size = 0.9) +
    geom_gene_arrow(data = sd_data,
                    arrowhead_height = unit(0, "mm"),
                    arrowhead_width = unit(0, "mm"),
                    arrow_body_height = unit(2.1, 'mm'),
                    aes(fill = color)) +
    scale_fill_manual(values = color_mapping) +
    geom_gene_arrow(fill = NA,
                    arrowhead_height = unit(3, "mm"), 
                    arrowhead_width = unit(1, "mm"),
                    arrow_body_height = unit(3, 'mm'),
                    size = 0.9) +
    facet_wrap(~ genome, scales = "free", ncol = 1) +
    labs(title = product) +
    theme_genes()
}
