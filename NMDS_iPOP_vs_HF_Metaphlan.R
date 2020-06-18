# 6/17/20
# Create NMDS graph for comparison ipop (non-hf) vs hf, ipop vs hf and ir vs is,
#   and ipop vs hf with nyha II-IV
# Filters out excluded samples (no metadata)

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/iPOP_vs_HF"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"
set.seed(10)
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

nmds_batch_unbatch_comp_hf <- function(met_lev, met_reg, met_lab, sing_lab, 
                                       nyhaExclude = c(), irisComp = F, desc, file_name) {
  print(met_lab)
  
  load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
  met_df <- met_taxa_df
  
  graph_dir2 <- file.path(graph_dir, paste(file_name, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # filter metadata
  met_metadata2 <- met_metadata %>% 
    filter(analysis != "excluded" &
             !(nyha %in% nyhaExclude))
  included_ids <- met_metadata2$study_id
  hf <- met_metadata2$hf
  names(hf) <- met_metadata2$study_id
  
  # setup ir/is
  if (irisComp) {
    met_metadata2 <- met_metadata2 %>% arrange(homa_ir)
    mid <- round(ncol(met_metadata2) / 2)
    ir_is <- c(rep("IS", mid), rep("IR", ncol(met_metadata2) - mid))
    names(ir_is) <- met_metadata2$study_id
  }

  # get taxa
  met_taxa_df <- met_df %>% 
    filter(batch == T &
             level == met_lev &
             study_id %in% included_ids) %>%
    mutate(hf = hf[study_id],
           taxa = str_replace(taxa, met_reg, ""))
  if (irisComp) met_taxa_df <- met_taxa_df %>%
    mutate(ir_is = ir_is[study_id])
  
  met_taxa_df <- met_taxa_df %>%
    # filter(taxa %in% top_N_phy) %>%
    spread(taxa, val)
  met_data <- met_taxa_df %>%
    select(-c(participant_id, study_id, collection_type, omni_comparison, batch, level, hf))
  if (irisComp) met_data <- met_data %>% select(-ir_is)
  met_data <- met_data %>% as.matrix()
  if (irisComp) {
    met_metadata <- met_taxa_df %>%
      select(ir_is, hf)
  } else {
    met_metadata <- met_taxa_df %>%
      select(hf)
  }
  
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
  
  if (irisComp) {
    gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = hf, color = hf, shape = ir_is)) + 
      geom_point(size = 3) +  # Set color to vary based on subject
      labs(title = paste(desc, met_lab, sep = " "))
  } else {
    gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = hf, color = hf)) + 
      geom_point(size = 3) +  # Set color to vary based on subject
      labs(title = paste(desc, met_lab, sep = " "))
  }
  plot(gg)
  ggsave(file.path(graph_dir2, paste(met_lab, "NMDS.pdf", sep = "_")), plot = last_plot())
}

# for each taxum
run_3 <- function(met_lab, sing_lab, met_lev, met_reg) {
  nmds_batch_unbatch_comp_hf(met_lev = met_lev, met_reg = met_reg, met_lab = met_lab, sing_lab = sing_lab,
                             desc = "iPOP vs Heart Failure Samples\nNMDS (Bray-Curtis Distance), ",
                             file_name = "iPOP_vs_Heart_Failure") # hf vs ipop
  nmds_batch_unbatch_comp_hf(met_lev = met_lev, met_reg = met_reg, met_lab = met_lab, sing_lab = sing_lab, irisComp = T,
                             desc = "iPOP vs Heart Failure Samples and IR vs IS\nNMDS (Bray-Curtis Distance), ",
                             file_name = "iPOP_vs_Heart_Failure_IR_vs_IS") # hf vs ipop, ir vs is
  nmds_batch_unbatch_comp_hf(met_lev = met_lev, met_reg = met_reg, met_lab = met_lab, sing_lab = sing_lab, nyhaExclude = 1,
                             desc = "iPOP vs Heart Failure Samples for NYHA 2-4\nNMDS (Bray-Curtis Distance), ",
                             file_name = "iPOP_vs_Heart_Failure_NYHA") # hf vs ipop, nyha 2-4 only
}

run_3("Phyla", "Phylum", "phylum", ".*p__")
run_3("Genera", "Genus", "genus", ".*g__")