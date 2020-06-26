# 6/25/20
# Create NMDS graph for Danish datasets, omni vs regular

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data/Tidy_Danish"
graph_dir <- "Graphs/Danish"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
set.seed(10)
# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

nmds_danish <- function(ds1) {
  print(ds1)
  
  graph_dir2 <- file.path(graph_dir, ds1)
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  desc <- paste(ds1, "NMDS, (Bray-Curtis Distance)")
  
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))

  # subset/setup data and metadata
  met_phyla_df <- tidy_df %>%
    spread(taxa, val)
  met_data <- met_phyla_df %>%
    arrange(comb_id) %>% # align data and metadata
    select(-comb_id) %>%
    as.matrix()
  met_metadata <- tidy_metadata %>%
    arrange(comb_id) # align data and metadata

  # get nmds
  print("Elbow")
  pdf(file.path(graph_dir2, paste(ds1, "Elbow_Plot.pdf", sep = "_")), width = 18, height = 12)
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
  
  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = TRIMESTER, color = TRIMESTER)) + 
    geom_point(size = 3) +
    labs(title = desc)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS.pdf", sep = "_")), plot = last_plot())
}

### comparing for metaphlan, phyla
nmds_danish("Metaphlan_Genus")
nmds_danish("uBiome")
nmds_danish("PCL")
nmds_danish("Metaphlan")
nmds_danish("Meta_Pathcoverage")