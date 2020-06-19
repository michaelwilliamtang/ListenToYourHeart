# 6/11/20
# Create NMDS graph for unbatched pcl, omni vs regular

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/Unbatched_Omni_vs_Regular"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"
set.seed(10)
desc <- "Omnigene vs Regular for Unbatched Samples\nNMDS (Bray-Curtis Distance)"
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

load(file.path(data_dir, paste("Tidy_Filtered_", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

nmds_omni_reg_comp <- function(pcl_lab, sing_lab, N, is_batch) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # get top pathways
  pcl_pathway_df <- pcl_df %>% 
    filter(omni_comparison == "comparison" &
             batch == is_batch)
  pcl_pathway_df$batched <- "unbatched"
  pcl_pathway_df$batched[which(pcl_pathway_df$batch)] <- "batched"
  phy <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$pathway[1:N]
  
  pcl_pathway_df <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_")) %>%
    spread(pathway, val)
  pcl_data <- pcl_pathway_df %>%
    select(top_N_phy) %>%
    as.matrix()
  pcl_metadata <- pcl_pathway_df %>%
    select(participant_id, collection_type)
  
  # get nmds
  print("Elbow")
  pdf(file.path(graph_dir2, paste(pcl_lab, "Elbow_Plot.pdf", sep = "_")), width = 18, height = 12)
  dimcheckMDS(pcl_data, distance = "bray", k = 10, trymax = 20,
              autotransform = TRUE)
  dev.off()
  print("NMDS")
  nmds <- metaMDS(pcl_data, distance = "bray", k = 2, trymax = 500, na.rm = T, autotransform = F)
  
  # envfit scores, write dimensions and vector wts
  nmds_scores <- as.data.frame(scores(nmds))  
  nmds_scores <- cbind(nmds_scores, pcl_metadata)
  write.table(nmds_scores, row.names = F, file = file.path(graph_dir2, "dimensions.txt"), sep = "\t", quote = FALSE)
  envf <- envfit(nmds, pcl_data, perm = 999, na.rm = T)
  vector.weights <-  scores(envf, display = "vectors")
  padjusted <- p.adjust(envf$vectors$pvals, method = "fdr")
  vector.weights <- cbind(analyte = rownames(vector.weights), vector.weights, padjust = padjusted)
  write.table(vector.weights, row.names = F, file = file.path(graph_dir2, "vector_weights.txt"), sep = "\t", quote = FALSE)
  
  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = participant_id, color = participant_id, shape = collection_type)) + 
    geom_point(size = 3) +  # Set color to vary based on subject
    labs(title = paste(desc, pcl_lab, sep = " "))
  plot(gg)
  ggsave(file.path(graph_dir2, paste(pcl_lab, "NMDS.pdf", sep = "_")), plot = last_plot())
}

### comparing for pcl
nmds_omni_reg_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, is_batch = F)

### same, but batched
graph_dir <- "Graphs/PCL/Batched_Omni_vs_Regular"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

nmds_omni_reg_comp(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, is_batch = T)