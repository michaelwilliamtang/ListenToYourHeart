# 6/11/20
# Compares batch omni vs unbatch omni vs batch reg vs omni reg PCL

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/Batched_vs_Unbatched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

batch_vs_unbatch_comp <- function(pcl_lab, sing_lab, N, comp_N) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # desc1 <- paste("Omnigene vs Regular and Batch Corrected vs Uncorrected\nTop",
  #   comp_N, pcl_lab, "By Collection Difference", sep = " ")
  desc2 <- paste("Omnigene vs Regular and Batch Corrected vs Uncorrected\nTop",
                 N, pcl_lab, "By Relative Abundance", sep = " ")
  
  pcl_pathway_df <- pcl_df %>% 
    filter(omni_comparison == "comparison")
  pcl_pathway_df$batched <- "unbatched"
  pcl_pathway_df$batched[which(pcl_pathway_df$batch)] <- "batched"
  phy <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$pathway[1:N]
  
  # plot top phyla
  pcl_pathway_df2 <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_")) %>%
    mutate(study_id = paste(study_id, batched, sep = "_")) %>%
    group_by(pathway, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(pcl_pathway_df)
  fp <- file.path(graph_dir2, paste("Top", N, pcl_lab, "Stack.pdf", sep = "_"))
  pdf(fp, width = 36, height = 24)
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
  # pcl_pathway_df2 <- pcl_pathway_df %>%
  #   filter(pathway %in% top_N_phy) %>%
  #   group_by(batched, pathway) %>%
  #   summarize(mean_val = mean(val)) %>%
  #   arrange(mean_val)
  # pdf(file.path(graph_dir2, paste(pcl_lab, "Aggregated.pdf", sep = "_")), width = 18, height = 12)
  # print(pcl_pathway_df2 %>% ggplot(aes(x = batched, y = mean_val, fill = pathway)) +
  #         geom_bar(position = "stack", stat = "identity") +
  #         xlab("Collection Type") +
  #         ylab("Relative Abundance") +
  #         theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
  #         guides(fill = guide_legend(title = sing_lab)))
  # dev.off()
}

### comparing for pcl, pathways
batch_vs_unbatch_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, comp_N = 20)