library(tidyverse)
data_dir <- "Data"
graph_dir <- "Graphs"
summarize <- dplyr::summarize

# read pcl
pcl_df <- read.table(file.path(data_dir, "hf-abundance-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
pcl_df <- abundance_df %>% as.matrix() %>% t() %>% as_tibble()
# View(pcl_df[,1:50])