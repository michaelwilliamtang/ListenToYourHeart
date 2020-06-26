# 6/25/20
# Compares omni vs reg dissimilarity for batch/unbatch PCL

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/Unbatched_Omni_vs_Regular_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"

load(file.path(data_dir, paste("Tidy_Filtered_", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

omni_reg_comp_dis <- function(pcl_lab, sing_lab, N, is_batch) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  pcl_pathway_df <- pcl_df %>% 
    filter(omni_comparison == "comparison" &
           batch == is_batch) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_"))
  path <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_path <- path$pathway[1:N]
  
  # matrix for unfiltered
  pcl_pathway_mat <- pcl_pathway_df %>%
    select(study_id, pathway, val) %>%
    spread(pathway, val)
  rownames(pcl_pathway_mat) <- pcl_pathway_mat$study_id
  pcl_pathway_mat <- pcl_pathway_mat %>% 
    select(-study_id) %>%
    as.matrix()
  dist <- pcl_pathway_mat %>% dist() %>% as.matrix()
  write.table(as.data.frame(dist), file.path(graph_dir2, paste(pcl_lab, "Dissimilarity_Matrix.tsv", sep = "_")),
              sep = "\t")
  
  # matrix for filtered
  pcl_pathway_mat <- pcl_pathway_df %>%
    filter(pathway %in% top_N_path) %>%
    select(study_id, pathway, val) %>%
    spread(pathway, val)
  rownames(pcl_pathway_mat) <- pcl_pathway_mat$study_id
  pcl_pathway_mat <- pcl_pathway_mat %>% 
    select(-study_id) %>%
    as.matrix()
  dist <- pcl_pathway_mat %>% dist() %>% as.matrix()
  write.table(as.data.frame(dist), file.path(graph_dir2, paste(pcl_lab, "Dissimilarity_Matrix_Filtered.tsv", sep = "_")),
            sep = "\t")
}

### comparing for pcl, pathway
omni_reg_comp_dis(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, is_batch = F)

### same, but batched
graph_dir <- "Graphs/PCL/Batched_Omni_vs_Regular_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for pcl, pathway
omni_reg_comp_dis(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, is_batch = T)