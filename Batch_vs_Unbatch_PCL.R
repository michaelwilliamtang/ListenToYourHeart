# 6/11/20
# Compares batch vs unbatch for omni or reg PCL

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/Omni_Batched_vs_Unbatched_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

batch_vs_unbatch_comp <- function(pcl_lab, sing_lab, N, comp_N, coll_type) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  pcl_pathway_df <- pcl_df %>% 
    filter(omni_comparison == "comparison" &
           collection_type == coll_type)
  pcl_pathway_df$batched <- "unbatched"
  pcl_pathway_df$batched[which(pcl_pathway_df$batch)] <- "batched"
  phy <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$pathway[1:N]
  
  # plot most diff phyla, omni vs reg
  comp_phyla <- pcl_pathway_df %>%
    group_by(pathway, batched) %>%
    summarize(mean_val = mean(val)) %>%
    spread(batched, mean_val) %>%
    mutate(collection_diff = batched - unbatched,
           collection_diff_abs = abs(collection_diff)) %>%
    arrange(desc(collection_diff_abs))
  diff_N_phy <- comp_phyla$pathway[1:comp_N]
  comp_phyla <- comp_phyla %>% filter(pathway %in% diff_N_phy)
  pdf(file.path(graph_dir2, paste(pcl_lab, "By_Collection_Diff.pdf", sep = "_")), width = 18, height = 12)
  print(comp_phyla %>% ggplot(aes(x = reorder(pathway, -collection_diff), y = collection_diff)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab(sing_lab) +
          ylab("Relative Abundance Diff (batched - unbatched)") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)))
  dev.off()
  
  # plot top phyla
  pcl_pathway_df2 <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    group_by(pathway, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(pcl_pathway_df)
  pdf(file.path(graph_dir2, paste("Top", N, pcl_lab, "Stack.pdf", sep = "_")), width = 18, height = 12)
  print(pcl_pathway_df2 %>% ggplot(aes(x = study_id, y = mean_val, fill = pathway)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Sample") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)))
  dev.off()
  pcl_pathway_df2 <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    group_by(batched, pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  pdf(file.path(graph_dir2, paste(pcl_lab, "Aggregated.pdf", sep = "_")), width = 18, height = 12)
  print(pcl_pathway_df2 %>% ggplot(aes(x = batched, y = mean_val, fill = pathway)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Collection Type") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)))
  dev.off()
}

### comparing for pcl, pathways
batch_vs_unbatch_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, comp_N = 20, coll_type = "omni")

### now, same but with regular
graph_dir <- "Graphs/PCL/Regular_Batched_vs_Unbatched"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for pcl, pathways
batch_vs_unbatch_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, comp_N = 20, coll_type = "regular")