# 6/25/20
# Create NMDS graphs for Danish datasets, colored by each metadatum,
#   with factor overlays

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
# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")
# library(ggvegan)

nmds_danish <- function(ds1, env_var, graph_dir_ds, quant, clean_name, env_vars) {
  print(env_var)
  
  graph_dir2 <- file.path(graph_dir_ds, env_var)
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  desc <- paste(ds1, "NMDS (Bray-Curtis Distance)")
  
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  # get rid of rows with missing metadata
  tidy_metadata[tidy_metadata == ""] <- NA
  missing <- tidy_metadata$comb_id[which(is.na(tidy_metadata[, env_var]))]
  tidy_df <- tidy_df %>% filter(!(comb_id %in% missing))
  tidy_metadata <- tidy_metadata %>% filter(!(comb_id %in% missing))
  
  # subset/setup data and metadata
  tidy_df <- tidy_df %>%
    spread(analyte, val)
  tidy_data <- tidy_df %>%
    arrange(comb_id) %>% # align data and metadata
    select(-comb_id) %>%
    as.matrix()
  tidy_metadata <- tidy_metadata %>%
    arrange(comb_id) # align data and metadata
  
  # rename
  if (quant) {
    tidy_metadata[, env_var] <- quantcut(as.numeric(tidy_metadata[, env_var]), q = 3, na.rm = T)
    levels(tidy_metadata[, env_var]) <- c("First Tertile", "Second Tertile", "Third Tertile")
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
  scrs <- as.data.frame(scores(nmds, display = "sites"))
  scrs <- cbind(scrs, tidy_metadata)
  tidy_meta <- as.data.frame(tidy_metadata[, env_var])
  # colnames(tidy_meta) <- paste0(clean_name, ", ") # include var prefix
  colnames(tidy_meta) <- "" # level only
  envf <- envfit(nmds, tidy_meta, choices = 1:3, perm = 999, na.rm = T)
  scrs2 <- as.data.frame(scores(envf, display = "factors"))
  scrs2 <- cbind(scrs2, lev = rownames(scrs2))
  
  # write scores with padjust
  write.table(scrs, row.names = F, file = file.path(graph_dir2, "dimensions.txt"), sep = "\t", quote = FALSE)
  padjusted <- p.adjust(envf$factors$pvals, method = "fdr")
  scrs3 <- cbind(scrs2, padjust = padjusted)
  write.table(scrs3, row.names = F, file = file.path(graph_dir2, "factor_weights.txt"), 
                          sep = "\t", quote = FALSE)
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS1, y = NMDS2, fill = !!sym(env_var), color = !!sym(env_var)), size = 3) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                     arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS1, y = NMDS2, label = lev), size = 3) +
    labs(title = desc,
         fill = clean_name,
         color = clean_name)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-2.pdf", sep = "_")), plot = last_plot())
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS2, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var)), size = 3) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS2, y = 0, yend = NMDS3),
                 arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS2, y = NMDS3, label = lev), size = 3) +
    labs(title = desc,
         fill = clean_name,
         color = clean_name)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS2-3.pdf", sep = "_")), plot = last_plot())
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS1, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var)), size = 3) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3),
                 arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS1, y = NMDS3, label = lev), size = 3) +
    labs(title = desc,
         fill = clean_name,
         color = clean_name)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-3.pdf", sep = "_")), plot = last_plot())
  
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
  clean_names <- scan(file.path("Metadata", "Danish_Env", paste0("Clean_Names_", ds1, ".txt")),
                      character(), quote = '', sep = "\t", quiet = T)
  names(clean_names) <- env_vars
  graph_dir_ds <- file.path(graph_dir, ds1)
  if (!dir.exists(graph_dir_ds)) dir.create(graph_dir_ds)
  for (env_var in env_vars) nmds_danish(ds1, env_var, graph_dir_ds, env_var %in% quant_vars,
                                        clean_names[env_var], env_vars)
}

### comparing for metaphlan, phyla
datasets <- c("Metaphlan_Genus", "uBiome", "PCL", "Metaphlan", "Meta_Pathcoverage")
for (ds in datasets) nmds_danish_all(ds)
