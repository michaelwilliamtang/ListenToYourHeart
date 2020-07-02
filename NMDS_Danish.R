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
library(gtools)

nmds_danish <- function(ds1, env_var, graph_dir_ds, quant) {
  print(env_var)
  
  graph_dir2 <- file.path(graph_dir_ds, env_var)
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  desc <- paste(ds1, "NMDS, (Bray-Curtis Distance)")
  
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))

  # subset/setup data and metadata
  met_phyla_df <- tidy_df %>%
    spread(analyte, val)
  tidy_data <- met_phyla_df %>%
    arrange(comb_id) %>% # align data and metadata
    select(-comb_id) %>%
    as.matrix()
  tidy_metadata <- tidy_metadata %>%
    arrange(comb_id) # align data and metadata
  tidy_metadata[tidy_metadata == ""] <- NA
  if (quant) {
    tidy_metadata[, env_var] <- quantcut(as.numeric(tidy_metadata[, env_var]), q = 3, na.rm = T)
  }

  # get nmds
  print("Elbow")
  pdf(file.path(graph_dir2, paste(ds1, "Elbow_Plot.pdf", sep = "_")), width = 18, height = 12)
  dimcheckMDS(tidy_data, distance = "bray", k = 10, trymax = 20,
              autotransform = TRUE)
  dev.off()
  print("NMDS")
  nmds <- metaMDS(tidy_data, distance = "bray", k = 3, trymax = 500, na.rm = T, autotransform = F)

  # envfit scores, write dimensions and vector wts
  nmds_scores <- as.data.frame(scores(nmds))
  nmds_scores <- cbind(nmds_scores, tidy_metadata)
  write.table(nmds_scores, row.names = F, file = file.path(graph_dir2, "dimensions.txt"), sep = "\t", quote = FALSE)
  # envf <- envfit(nmds, tidy_metadata %>% select(env_vars), perm = 999, na.rm = T)
  # vector.weights <-  scores(envf, display = "vectors")
  # padjusted <- p.adjust(envf$vectors$pvals, method = "fdr")
  # vector.weights <- cbind(analyte = rownames(vector.weights), vector.weights, padjust = padjusted)
  # write.table(vector.weights, row.names = F, file = file.path(graph_dir2, "vector_weights.txt"), sep = "\t",
  # quote = FALSE)

  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = !!sym(env_var), color = !!sym(env_var))) +
    geom_point(size = 3) +
    labs(title = desc)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-2.pdf", sep = "_")), plot = last_plot())

  gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var))) +
    geom_point(size = 3) +
    labs(title = desc)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-3.pdf", sep = "_")), plot = last_plot())

  gg <- ggplot(nmds_scores, aes(x = NMDS2, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var))) +
    geom_point(size = 3) +
    labs(title = desc)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS2-3.pdf", sep = "_")), plot = last_plot())

  
  ### for old ubiome, no overlay:
  # gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  #   geom_point(size = 3) +
  #   labs(title = desc)
  # plot(gg)
  # ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-2.pdf", sep = "_")), plot = last_plot())
  # 
  # gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS3)) +
  #   geom_point(size = 3) +
  #   labs(title = desc)
  # plot(gg)
  # ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-3.pdf", sep = "_")), plot = last_plot())
  # 
  # gg <- ggplot(nmds_scores, aes(x = NMDS2, y = NMDS3)) +
  #   geom_point(size = 3) +
  #   labs(title = desc)
  # plot(gg)
  # ggsave(file.path(graph_dir2, paste(ds1, "NMDS2-3.pdf", sep = "_")), plot = last_plot())
}

nmds_danish_all <- function(ds1) {
  print(ds1)
  env_vars <- scan(file.path("Metadata", "Danish_Env", paste0("Env_", ds1, ".txt")),
                   character(), quote = '', sep = "\t", quiet = T) %>%
              str_replace_all("\\(|\\)", "\\.")
  quant_vars <- scan(file.path("Metadata", "Danish_Env", paste0("Quant_", ds1, ".txt")),
                   character(), quote = '', sep = "\t", quiet = T) %>%
              str_replace_all("\\(|\\)", "\\.")
  graph_dir_ds <- file.path(graph_dir, ds1)
  if (!dir.exists(graph_dir_ds)) dir.create(graph_dir_ds)
  for (env_var in env_vars) nmds_danish(ds1, env_var, graph_dir_ds, env_var %in% quant_vars)
}

### comparing for metaphlan, phyla
datasets <- c("Metaphlan_Genus", "uBiome", "PCL", "Metaphlan", "Meta_Pathcoverage")
for (ds in datasets) nmds_danish_all(ds)
