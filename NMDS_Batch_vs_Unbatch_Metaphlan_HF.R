# 6/12/20
# Create NMDS graph for comparison of batch vs unbatch for hf patients in metaphlan

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/HF_Batched_vs_Unbatched_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"
set.seed(10)
desc <- "Batch Corrected vs Uncorrected, Heart Failure Samples\nNMDS (Bray-Curtis Distance)"
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
metaphlan_df <- met_taxa_df

nmds_batch_unbatch_comp_hf <- function(met_lev, met_reg, met_lab, sing_lab, N, hf_tok) {
  print(met_lab)
  
  graph_dir2 <- file.path(graph_dir, paste("Top", N, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # setup hf check
  is_hf <- met_metadata$hf
  names(is_hf) <- met_metadata$study_id
  
  # get top taxa
  met_phyla_df <- metaphlan_df %>% 
    filter(level == met_lev &
           is_hf[study_id] == hf_tok) %>%
    mutate(taxa = str_replace(taxa, met_reg, ""))
  met_phyla_df$batched <- "unbatched"
  met_phyla_df$batched[which(met_phyla_df$batch)] <- "batched"
  phy <- met_phyla_df %>%
    group_by(taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$taxa[1:N]
  
  # subset/setup data and metadata
  met_phyla_df <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, batched, sep = "_")) %>%
    spread(taxa, val)
  met_data <- met_phyla_df %>%
    select(top_N_phy) %>%
    as.matrix()
  met_metadata <- met_phyla_df %>%
    select(participant_id, batched)

  # get nmds
  print("Elbow")
  pdf(file.path(graph_dir2, paste(met_lab, "Elbow_Plot.pdf", sep = "_")), width = 18, height = 12)
  dimcheckMDS(met_data, distance = "bray", k = 10, trymax = 20,
              autotransform = TRUE)
  dev.off()
  print("NMDS")
  nmds <- metaMDS(met_data, distance = "bray", k = 2, trymax = 500, na.rm = T, autotransform = F)
  
  # envfit scores, write dimensions and vector wts
  nmds_scores <- as.data.frame(scores(nmds))  
  nmds_scores <- cbind(nmds_scores, met_metadata)
  write.table(nmds_scores, row.names = F, file = file.path(graph_dir2, "dimensions.txt"), sep = "\t", quote = FALSE)
  envf <- envfit(nmds, met_data, perm = 999, na.rm = T)
  vector.weights <-  scores(envf, display = "vectors")
  padjusted <- p.adjust(envf$vectors$pvals, method = "fdr")
  vector.weights <- cbind(analyte = rownames(vector.weights), vector.weights, padjust = padjusted)
  write.table(vector.weights, row.names = F, file = file.path(graph_dir2, "vector_weights.txt"), sep = "\t", quote = FALSE)
  
  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = participant_id, color = participant_id, shape = batched)) + 
    geom_point(size = 3) +  # Set color to vary based on subject
    labs(title = paste(desc, sing_lab, "Level", sep = " "))
  plot(gg)
  ggsave(file.path(graph_dir2, paste(met_lab, "NMDS.pdf", sep = "_")), plot = last_plot())
}

### comparing for metaphlan, phyla
nmds_batch_unbatch_comp_hf(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, hf_tok = "1")

### comparing for metaphlan, genera
nmds_batch_unbatch_comp_hf(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, hf_tok = "1")