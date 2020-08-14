# 6/25/20
# Create NMDS graphs for Danish datasets, colored by each metadatum,
#   with factor overlays, with an option to aggregate per participant

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- file.path("Data", "Tidy_Danish")
graph_dir <- file.path("Graphs", "Danish")
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

nmds_danish <- function(ds1, env_var, graph_dir_ds, quant, clean_name, env_vars, aggregated) {
  print(env_var)
  
  graph_dir2 <- file.path(graph_dir_ds, env_var)
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  desc <- paste(ds1, "NMDS (Bray-Curtis Distance)")
  if (aggregated) desc <- paste0(desc, ", Aggregated Samples")
  
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  # get rid of rows with missing metadata
  tidy_metadata[tidy_metadata == ""] <- NA
  missing <- tidy_metadata$comb_id[which(is.na(tidy_metadata[, env_var]))]
  tidy_df <- tidy_df %>% filter(!(comb_id %in% missing))
  tidy_metadata <- tidy_metadata %>% filter(!(comb_id %in% missing))
  
  # aggregate per participant if specified
  if (aggregated) {
    ids <- tidy_metadata$PARTICIPANT_ID
    names(ids) <- tidy_metadata$comb_id
    tidy_df <- tidy_df %>%
      mutate(comb_id = ids[comb_id]) %>% # use part id as comb_id, so only unique per participant
      group_by(comb_id, analyte) %>%
      summarize(val = mean(val), .groups = "drop")
    tidy_metadata <- tidy_metadata %>%
      mutate(comb_id = PARTICIPANT_ID) %>%
      select(!! sym(env_var), comb_id, PARTICIPANT_ID)
    tidy_metadata[, env_var] <- as.numeric(tidy_metadata[, env_var]) # not reliably num
    tidy_metadata <- aggregate(tidy_metadata[, env_var], list(comb_id = tidy_metadata$comb_id), mean)
    
    # minor post-aggregate cleanup
    colnames(tidy_metadata)[2] = env_var # analyte is renamed as "x" by aggregate
    tidy_metadata$PARTICIPANT_ID <- tidy_metadata$comb_id # other cols get removed by aggregate
  }
  
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
  if (quant || env_var == "TRIMESTER") { # trimester data changes over visits, after aggregate may not be integer
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
  
  sz <- 4
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS1, y = NMDS2, fill = !!sym(env_var), color = !!sym(env_var)), size = sz) +
    geom_text(aes(x = NMDS1, y = NMDS2, label = PARTICIPANT_ID), size = sz) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
                     arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS1, y = NMDS2, label = lev), size = sz) +
    labs(title = desc,
         fill = clean_name,
         color = clean_name)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS1-2.pdf", sep = "_")), plot = last_plot())
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS2, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var)), size = sz) +
    geom_text(aes(x = NMDS2, y = NMDS3, label = PARTICIPANT_ID), size = sz) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS2, y = 0, yend = NMDS3),
                 arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS2, y = NMDS3, label = lev), size = sz) +
    labs(title = desc,
         fill = clean_name,
         color = clean_name)
  plot(gg)
  ggsave(file.path(graph_dir2, paste(ds1, "NMDS2-3.pdf", sep = "_")), plot = last_plot())
  
  gg <- ggplot(scrs) +
    geom_point(aes(x = NMDS1, y = NMDS3, fill = !!sym(env_var), color = !!sym(env_var)), size = sz) +
    geom_text(aes(x = NMDS1, y = NMDS3, label = PARTICIPANT_ID), size = sz) +
    coord_fixed() + # fix aspect ratio
    geom_segment(data = scrs2,
                 aes(x = 0, xend = NMDS1, y = 0, yend = NMDS3),
                 arrow = arrow(length = unit(0.25, "cm")), color = "darkgray") +
    geom_text(data = scrs2, aes(x = NMDS1, y = NMDS3, label = lev), size = sz) +
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

nmds_danish_all <- function(ds1, aggregated = F) {
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
  dir_name <- ds1
  if (aggregated) dir_name <- paste0(dir_name, "_Aggregated")
  graph_dir_ds <- file.path(graph_dir, dir_name)
  if (!dir.exists(graph_dir_ds)) dir.create(graph_dir_ds)
  for (env_var in env_vars) {
    tryCatch({
      nmds_danish(ds1, env_var, graph_dir_ds, env_var %in% quant_vars,
                  clean_names[env_var], env_vars, aggregated)
    },
    error = function(cond) {
      message(paste0("Error with ", env_var, " in dataset ", ds1))
      message(cond)
    }
    )
  }
}

### per datset
datasets <- c("Metaphlan_Genus", "uBiome", "PCL", "Metaphlan", "Meta_Pathcoverage")
for (ds in datasets) nmds_danish_all(ds, T)
# for (ds in datasets) nmds_danish_all(ds, F)
