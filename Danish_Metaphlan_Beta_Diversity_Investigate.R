# 8/18/20
# Calculates beta diversity for metaphlan at genus level and plots by environmental metadatum

library(tidyverse)
library(vegan)
library(gtools)
data_dir <- file.path("Data", "Tidy_Danish")
graph_dir <- file.path("Graphs", "Danish", "Metaphlan_Beta_Diversity")
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize

beta_inv_all <- function() {
  # load
  ds1 <- "Metaphlan_Genus"
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  # spread into matrix
  tidy_df <- tidy_df %>% 
    spread(analyte, val)
  tidy_mat <- tidy_df %>%
    select(-comb_id) %>%
    as.matrix()
  
  # partition env metadata into tertiles
  env_vars <- scan(file.path("Metadata", "Danish_Env", paste0("Env_", ds1, ".txt")),
                   character(), quote = '', sep = "\t", quiet = T) %>%
    str_replace_all("\\(|\\)", "\\.")
  quant_vars <- scan(file.path("Metadata", "Danish_Env", paste0("Quant_", ds1, ".txt")),
                     character(), quote = '', sep = "\t", quiet = T) %>%
    str_replace_all("\\(|\\)", "\\.")
  clean_names <- scan(file.path("Metadata", "Danish_Env", paste0("Clean_Names_", ds1, ".txt")),
                      character(), quote = '', sep = "\t", quiet = T)
  names(clean_names) <- env_vars
  for (env_var in env_vars) {
    if (env_var %in% quant_vars) tidy_metadata[, env_var] <- quantcut(as.numeric(tidy_metadata[, env_var]), q = 3, na.rm = T)
  }
  tidy_meta <- tidy_metadata %>% 
    arrange(comb_id) # align with data
  
  # clean spaces to NA (betadisper does not like spaces)
  tidy_meta[tidy_meta == ""] = NA
  
  # remove NAs
  # rm_cols <- is.na(tidy_meta$TRIMESTER)
  # tidy_meta <- tidy_meta[!rm_cols,]
  # tidy_mat <- tidy_mat[!rm_cols,]
  
  # calc beta diversity indices
  beta_vec <- betadiver(tidy_mat, "z")
  for (env_var in env_vars) beta_inv_metadatum(beta_vec, tidy_meta, env_var, clean_names[env_var])
}

beta_inv_metadatum <- function(beta_vec, tidy_meta, env_var, clean_name) {
  print(env_var)
  
  # add metadata
  beta_plot <- with(tidy_meta, betadisper(beta_vec,  eval(parse(text = env_var)))) # treat string as variable name
  
  # subdir
  graph_dir2 <- file.path(graph_dir, env_var)
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # plot
  pdf(file.path(graph_dir2, paste0("Beta_Diversity_Boxplot_", env_var, ".pdf")))
  boxplot(beta_plot, main = paste0("Beta Diversity Distance by ", clean_name), xlab = clean_name)
  dev.off()
  pdf(file.path(graph_dir2, paste0("Beta_Diversity_Plot_", env_var, ".pdf")))
  plot(beta_plot, main = paste0("Beta Diversity Distance by ", clean_name), xlab = clean_name)
  dev.off()
}

beta_inv_all()
