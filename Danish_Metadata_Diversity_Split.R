# 8/18/20
# Plot analyte-metadatum pairs from Danish metaphlan genus diversity data

library(tidyverse)
data_dir <- file.path("Data", "Tidy_Danish")
graph_dir <- file.path("Graphs", "Danish")
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize

dci_diversity_metadatum <- function(data_df, metadata_df, cname, meta1, graph_dir3) {
  # convert metadata (default factor for some reason)
  metadata_df[,meta1] <- as.numeric(metadata_df[,meta1])
  
  all_analytes <- data_df$analyte %>% unique()

  # plot those top N analytes with selected metadatum
  for (sel in all_analytes) {
    tryCatch({
      # get metadatum-analyte pair
      sel_data <- data_df %>% filter(analyte == sel) %>% select(comb_id, val)
      sel_meta <- metadata_df %>% select(comb_id, !! sym(meta1), PARTICIPANT_ID)
      
      # filter missing metadata, merge
      valid_ids <- sel_meta$comb_id[!is.na(sel_meta[,meta1])]
      sel_data <- sel_data %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      sel_meta <- sel_meta %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      print(table(sel_data$comb_id == sel_meta$comb_id))
      comb_sel <- cbind(sel_data %>% select(val, comb_id), meta = sel_meta[, meta1], id = sel_meta$PARTICIPANT_ID)
      
      # test sum == 0
      comb_sum <- comb_sel$meta %>% sum()
      if (comb_sum == 0) {
        message(paste0("All zeroes with ", sel, " for ", cname))
        next
      }
      
      med_meta <- comb_sel$meta %>% median()
      comb_sel <- comb_sel %>% mutate(half = ifelse(meta <= med_meta, "low", "high"))
      comb_sel$half <- factor(comb_sel$half, ordered = T, levels = c("low", "high"))
      
      # plot
      desc <- paste(cname, "vs", sel, "Split Distribution", sep = " ")
      print(desc)
      gg <- ggplot(comb_sel, aes(x = half, y = val, fill = half)) +
        geom_violin(trim = F) +
        # geom_text(label = comb_sel$comb_id, 
        #           position = position_dodge(width = 1)) +
        geom_jitter(position = position_jitter(0.2), aes(color = id)) +
        xlab(cname) +
        ylab(sel) +
        labs(title = desc, color = "Participant", fill = "Half")
      
      plot(gg)
      ggsave(file.path(graph_dir3, paste(sel, meta1, "Distributions.pdf", sep = "_")), plot = last_plot())
    },
    error = function(cond) {
      message(paste0("Error with ", sel, " for ", cname))
      message(cond)
    },
    warning = function(cond) {
      message(paste0("Warning with ", sel, " for ", cname))
      message(cond)
    }
    )
  }
}


danish_diversity_all <- function(selected = NA) {
  ds1 <- "Metaphlan_Diversity"
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
  
  if (is.na(selected)) selected = env_vars 
  for (sel in selected) {
    tryCatch({
      load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
      graph_dir2 <- file.path(graph_dir, paste(ds1, "Split", sep = "_"))
      if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
      
      graph_dir3 <- file.path(graph_dir2, sel)
      if (!dir.exists(graph_dir3)) dir.create(graph_dir3)
      dci_diversity_metadatum(tidy_df, tidy_metadata, clean_names[sel], sel, graph_dir3)
    },
    error = function(cond) {
      message(paste0("Error with ", env_var, " in dataset ", ds1))
      message(cond)
    }
    )
  }
}

### per dataset
danish_diversity_all()

