# 6/17/20
# Create NMDS graph for comparison ipop (non-hf) vs hf, ipop vs hf and ir vs is,
#   and ipop vs hf with nyha II-IV
# Filters out excluded samples (no metadata)

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/PCL/iPOP_vs_HF_Scaled_Filtered"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "pcl"
set.seed(10)
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

nmds_batch_unbatch_comp_hf <- function(pcl_lab, sing_lab, nyhaExclude = c(), irisComp = F, desc, file_name) {
  load(file.path(data_dir, paste("Tidy_Scaled_Filtered_", ds1, ".RData", sep = "")))
  pcl_df <- pcl_pathway_df
  
  graph_dir2 <- file.path(graph_dir, paste(file_name, pcl_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # filter metadata
  pcl_metadata2 <- pcl_metadata %>% 
    filter(analysis != "excluded" &
             !(nyha %in% nyhaExclude))
  included_ids <- pcl_metadata2$study_id
  hf <- pcl_metadata2$hf
  names(hf) <- pcl_metadata2$study_id
  
  # setup ir/is
  if (irisComp) {
    pcl_metadata2 <- pcl_metadata2 %>% arrange(homa_ir)
    mid <- round(ncol(pcl_metadata2) / 2)
    ir_is <- c(rep("IS", mid), rep("IR", ncol(pcl_metadata2) - mid))
    names(ir_is) <- pcl_metadata2$study_id
  }

  # get pathways
  pcl_pathway_df <- pcl_df %>% 
    filter(batch == T & 
             study_id %in% included_ids) %>%
    mutate(hf = hf[study_id])
  if (irisComp) pcl_pathway_df <- pcl_pathway_df %>%
    mutate(ir_is = ir_is[study_id])
  
  pcl_pathway_df <- pcl_pathway_df %>%
    # filter(pathway %in% top_N_phy) %>%
    spread(pathway, val)
  pcl_data <- pcl_pathway_df %>%
    select(-c(participant_id, study_id, collection_type, omni_comparison, batch, hf))
  if (irisComp) pcl_data <- pcl_data %>% select(-ir_is)
  pcl_data <- pcl_data %>%
    as.matrix()
  if (irisComp) {
    pcl_metadata <- pcl_pathway_df %>%
      select(ir_is, hf)
  } else {
    pcl_metadata <- pcl_pathway_df %>%
      select(hf)
  }
  
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
  # envf <- envfit(nmds, pcl_data, perm = 999, na.rm = T)
  # vector.weights <-  scores(envf, display = "vectors")
  # padjusted <- p.adjust(envf$vectors$pvals, method = "fdr")
  # vector.weights <- cbind(analyte = rownames(vector.weights), vector.weights, padjust = padjusted)
  # write.table(vector.weights, row.names = F, file = file.path(graph_dir2, "vector_weights.txt"), sep = "\t", quote = FALSE)
  
  if (irisComp) {
    gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = hf, color = hf, shape = ir_is)) + 
      geom_point(size = 3) +  # Set color to vary based on subject
      labs(title = paste(desc, pcl_lab, sep = " "))
  } else {
    gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = hf, color = hf)) + 
      geom_point(size = 3) +  # Set color to vary based on subject
      labs(title = paste(desc, pcl_lab, sep = " "))
  }
  plot(gg)
  ggsave(file.path(graph_dir2, paste(pcl_lab, "NMDS.pdf", sep = "_")), plot = last_plot())
  
  # if (irisComp) {
  #   gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS3, fill = hf, color = hf, shape = ir_is)) + 
  #     geom_point(size = 3) +  # Set color to vary based on subject
  #     labs(title = paste(desc, pcl_lab, sep = " "))
  # } else {
  #   gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS3, fill = hf, color = hf)) + 
  #     geom_point(size = 3) +  # Set color to vary based on subject
  #     labs(title = paste(desc, pcl_lab, sep = " "))
  # }
  # plot(gg)
  # ggsave(file.path(graph_dir2, paste(pcl_lab, "NMDS1-3.pdf", sep = "_")), plot = last_plot())
  # 
  # if (irisComp) {
  #   gg <- ggplot(nmds_scores, aes(x = NMDS2, y = NMDS3, fill = hf, color = hf, shape = ir_is)) + 
  #     geom_point(size = 3) +  # Set color to vary based on subject
  #     labs(title = paste(desc, pcl_lab, sep = " "))
  # } else {
  #   gg <- ggplot(nmds_scores, aes(x = NMDS2, y = NMDS3, fill = hf, color = hf)) + 
  #     geom_point(size = 3) +  # Set color to vary based on subject
  #     labs(title = paste(desc, pcl_lab, sep = " "))
  # }
  # plot(gg)
  # ggsave(file.path(graph_dir2, paste(pcl_lab, "NMDS2-3.pdf", sep = "_")), plot = last_plot())
}

nmds_batch_unbatch_comp_hf(pcl_lab = "Pathways", sing_lab = "Pathway", 
                           desc = "iPOP vs Heart Failure Samples\nNMDS (Bray-Curtis Distance), ",
                           file_name = "iPOP_vs_Heart_Failure") # hf vs ipop
nmds_batch_unbatch_comp_hf(pcl_lab = "Pathways", sing_lab = "Pathway", irisComp = T,
                           desc = "iPOP vs Heart Failure Samples and IR vs IS\nNMDS (Bray-Curtis Distance), ",
                           file_name = "iPOP_vs_Heart_Failure_IR_vs_IS") # hf vs ipop, ir vs is
nmds_batch_unbatch_comp_hf(pcl_lab = "Pathways", sing_lab = "Pathway", nyhaExclude = 1,
                           desc = "iPOP vs Heart Failure Samples for NYHA 2-4\nNMDS (Bray-Curtis Distance), ",
                           file_name = "iPOP_vs_Heart_Failure_NYHA") # hf vs ipop, nyha 2-4 only