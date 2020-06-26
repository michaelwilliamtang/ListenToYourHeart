# 6/10/20
# Compares omni vs reg for batch PCL

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/Unbatched_Omni_vs_Regular_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"

load(file.path(data_dir, paste("Tidy_Filtered", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

omni_reg_comp <- function(pcl_lab, sing_lab, N, comp_N, is_batch) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  desc1 <- paste("Omnigene vs Regular, Top", comp_N, pcl_lab, "By Collection Difference", sep = " ")
  desc2 <- paste("Omnigene vs Regular, Top", N, pcl_lab, "By Relative Abundance", sep = " ")
  
  pcl_pathway_df <- pcl_df %>% 
    filter(omni_comparison == "comparison" &
           batch == is_batch)
  phy <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$pathway[1:N]
  
  # plot most diff phyla, omni vs reg
  comp_phyla <- pcl_pathway_df %>%
    group_by(pathway, collection_type) %>%
    summarize(mean_val = mean(val)) %>%
    spread(collection_type, mean_val) %>%
    mutate(collection_diff = omni - regular,
           collection_diff_abs = abs(collection_diff)) %>%
    arrange(desc(collection_diff_abs))
  diff_N_phy <- comp_phyla$pathway[1:comp_N]
  comp_phyla <- comp_phyla %>% filter(pathway %in% diff_N_phy)
  fp <- file.path(graph_dir2, paste(pcl_lab, "By_Collection_Diff.pdf", sep = "_"))
  pdf(fp,  width = 18, height = 12)
  gg <- comp_phyla %>% ggplot(aes(x = reorder(pathway, -collection_diff), y = collection_diff)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab(sing_lab) +
          ylab("Relative Abundance Diff (omni - reg)") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
         labs(title = desc1)
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
  
  # plot top phyla
  pcl_pathway_df2 <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_")) %>%
    group_by(pathway, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(pcl_pathway_df)
  fp <- file.path(graph_dir2, paste("Top", N, pcl_lab, "Stack.pdf", sep = "_"))
  pdf(fp, width = 18, height = 12)
  gg <- pcl_pathway_df2 %>% ggplot(aes(x = study_id, y = mean_val, fill = pathway)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Sample") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)) +
          labs(title = desc2)
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
  pcl_pathway_df2 <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    group_by(collection_type, pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  fp <- file.path(graph_dir2, paste(pcl_lab, "Aggregated.pdf", sep = "_"))
  pdf(fp,  width = 18, height = 12)
  gg <- pcl_pathway_df2 %>% ggplot(aes(x = collection_type, y = mean_val, fill = pathway)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Collection Type") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)) +
          labs(title = paste(desc2, ", Aggregated", sep = ""))
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
}

### comparing for pcl, pathways
omni_reg_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, comp_N = 20, is_batch = F)

### same, but batched
graph_dir <- "Graphs/PCL/Batched_Omni_vs_Regular_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for pcl, pathways
omni_reg_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, comp_N = 20, is_batch = T)