# 6/12/20
# Create NMDS graph for comparison of batch vs unbatch for hf patients in PCL

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/HF_Batched_vs_Unbatched_Scaled_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"
set.seed(10)
desc <- "Batch Corrected vs Uncorrected, Heart Failure Samples\nNMDS (Bray-Curtis Distance)"
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

load(file.path(data_dir, paste("Tidy_Scaled_Filtered_", ds1, ".RData", sep = "")))
pcl_df <- pcl_pathway_df

nmds_batch_unbatch_comp_hf <- function(pcl_lab, sing_lab, N, hf_tok) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # setup hf check
  is_hf <- pcl_metadata$hf
  names(is_hf) <- pcl_metadata$study_id
  
  # get top pathways
  pcl_pathway_df <- pcl_df %>% 
    filter(is_hf[study_id] == hf_tok)
  pcl_pathway_df$batched <- "unbatched"
  pcl_pathway_df$batched[which(pcl_pathway_df$batch)] <- "batched"
  phy <- pcl_pathway_df %>%
    group_by(pathway) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$pathway[1:N]
  
  pcl_pathway_df <- pcl_pathway_df %>%
    filter(pathway %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, batched, sep = "_")) %>%
    spread(pathway, val)
  pcl_data <- pcl_pathway_df %>%
    select(top_N_phy) %>%
    as.matrix()
  pcl_metadata <- pcl_pathway_df %>%
    select(participant_id, batched)
  
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
  
  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = participant_id, color = participant_id, shape = batched)) + 
    geom_point(size = 3) +  # Set color to vary based on subject
    labs(title = paste(desc, pcl_lab, sep = " "))
  plot(gg)
  ggsave(file.path(graph_dir2, paste(pcl_lab, "NMDS.pdf", sep = "_")), plot = last_plot())
}

### comparing for pcl
nmds_batch_unbatch_comp_hf(pcl_lab = "Pathways", sing_lab = "Pathway", N = 40, hf_tok = 1)