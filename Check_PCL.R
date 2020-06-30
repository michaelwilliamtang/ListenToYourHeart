# 6/11/20
# Checks summation of unlogged pcl values and graphs distributions

library(tidyverse)
data_dir <- "Data/Tidy"
summarize <- dplyr::summarize
ds1 <- "pcl"

check_distribution <- function(is_batch = F, coll_type = NA, file_prefix = "") {
  graph_dir2 <- file.path("Graphs", "PCL", "Check_Distribution")
  if (file_prefix != "") graph_dir2 <- paste(graph_dir2, file_prefix, sep = "_")
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  if (file_prefix != "") {
    load(file.path(data_dir, paste("Tidy_", file_prefix, "_", ds1, ".RData", sep = "")))
  } else {
    load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  }
  pcl_df <- pcl_pathway_df
  
  if (!is.na(is_batch)) pcl_df <- pcl_df %>%
    filter(batch == is_batch)
  if (!is.na(coll_type)) pcl_df <- pcl_df %>%
    filter(collection_type == coll_type)
  
  pcl_pathway_df <- pcl_df %>%
    filter(omni_comparison == "comparison")
  pcl_pathway_df$batched <- "unbatched"
  pcl_pathway_df$batched[which(pcl_pathway_df$batch)] <- "batched"
  
  # unlog
  # pcl_pathway_df$val <- 2 ^ pcl_pathway_df$val - 1
  
  check_df <- pcl_pathway_df %>% 
    group_by(batched, collection_type, study_id) %>%
    summarize(total = sum(val)) %>%
    arrange(desc(total))
  write_tsv(check_df, file.path(graph_dir2, paste("Totals.txt", sep = "_")))
  
  check_df <- pcl_pathway_df %>% 
    group_by(batched, collection_type, study_id, pathway) %>%
    summarize(total = sum(val)) %>%
    arrange(desc(total))
  write_tsv(check_df, file.path(graph_dir2, paste("Pathway_Totals.txt", sep = "_")))
  
  pdf(file.path(graph_dir2, paste("Check_Distribution.pdf", sep = "_")), width = 18, height = 12)
  print(pcl_pathway_df %>% ggplot(aes(x = val, fill = pathway)) +
    geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity") +
    theme(legend.position = "none"))
  dev.off()
}

check_distribution(file_prefix = "Scaled_Filtered", is_batch = T)
check_distribution(file_prefix = "Scaled", is_batch = T)
check_distribution(file_prefix = "Filtered", is_batch = T)
check_distribution(is_batch = T)