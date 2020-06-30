# 6/25/20
# Compares omni vs reg dissimilarity for batch/unbatch metaphlan

library(tidyverse)
library(ecodist)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/Unbatched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
metaphlan_df <- met_taxa_df

omni_reg_comp_dis <- function(met_lev, met_reg, met_lab, sing_lab, N, is_batch) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  met_phyla_df <- metaphlan_df %>% 
    filter(omni_comparison == "comparison" &
           level == met_lev &
           batch == is_batch) %>%
    mutate(taxa = str_replace(taxa, met_reg, ""),
           study_id = paste(study_id, collection_type, sep = "_"))
  phy <- met_phyla_df %>%
    group_by(taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$taxa[1:N]
  
  # matrix for unfiltered
  met_phyla_mat <- met_phyla_df %>%
    select(study_id, taxa, val) %>%
    spread(taxa, val)
  rownames(met_phyla_mat) <- met_phyla_mat$study_id
  met_phyla_mat <- met_phyla_mat %>% 
    select(-study_id) %>%
    as.matrix()
  dist <- met_phyla_mat %>% bcdist() %>% as.matrix()
  write.table(as.data.frame(dist), file.path(graph_dir2, paste(met_lab, "Dissimilarity_Matrix.tsv", sep = "_")),
              sep = "\t")
  
  # matrix for filtered
  met_phyla_mat <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    select(study_id, taxa, val) %>%
    spread(taxa, val)
  rownames(met_phyla_mat) <- met_phyla_mat$study_id
  met_phyla_mat <- met_phyla_mat %>% 
    select(-study_id) %>%
    as.matrix()
  dist <- met_phyla_mat %>% bcdist() %>% as.matrix()
  write.table(as.data.frame(dist), file.path(graph_dir2, paste(met_lab, "Dissimilarity_Matrix_Subset.tsv", sep = "_")),
            sep = "\t")
}

### comparing for metaphlan, phyla
omni_reg_comp_dis(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, is_batch = F)

### comparing for metaphlan, genera
omni_reg_comp_dis(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, is_batch = F)

### same, but batched
graph_dir <- "Graphs/Metaphlan/Batched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for metaphlan, phyla
omni_reg_comp_dis(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, is_batch = T)

### comparing for metaphlan, genera
omni_reg_comp_dis(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, is_batch = T)