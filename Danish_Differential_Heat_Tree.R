# 8/18/20
# Plots differential heat tree for metaphlan

library(tidyverse)
library(vegan)
library(gtools)
library(metacoder)
data_dir <- file.path("Data", "Tidy_Danish")
summarize <- dplyr::summarize
graph_dir <- file.path("Graphs", "Danish", "Metaphlan_Beta_Diversity_Tree")
if (!dir.exists(graph_dir)) dir.create(graph_dir)
set.seed(10)

diff_heat_trees <- function(ds1) {
  # load and clean
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  tidy_data <- tidy_df %>% 
    spread(analyte, val)
  cnames <- tidy_data$comb_id
  tidy_data <- tidy_data %>% select(-comb_id) %>% 
    t() %>% as.data.frame()
  colnames(tidy_data) <- cnames
  tidy_data$lineage <- rownames(tidy_data)
  tidy_data$lineage <- gsub("^unclassified", "k__unclassified", tidy_data$lineage, perl = T) # regex looks for level char
  tidy_metadata <- tidy_metadata %>% arrange(comb_id) # important to have same order when passing metadata as grouping
  
  # refactor NA and empty group names
  tidy_metadata[tidy_metadata == ""] = "UNKNOWN"
  tidy_metadata[tidy_metadata == "?"] = "UNKNOWN"
  tidy_metadata[is.na(tidy_metadata)] = "UNKNOWN"
  
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
    if (env_var %in% quant_vars) tidy_metadata[, env_var] <- quantcut(as.numeric(tidy_metadata[, env_var]), q = 2, na.rm = T)
  }
  tidy_metadata <- tidy_metadata %>% na.omit()
  
  # parse into taxmap
  obj <- parse_tax_data(tidy_data,
                        class_cols = "lineage",
                        class_sep = "|",
                        class_regex = "^(.+)__(.+)$",
                        class_key = c(tax_rank = "info",
                                      tax_name = "taxon_name")) # 2 capture groups
  obj$data$type_abund <- calc_taxon_abund(obj, "tax_data",
                                          cols = tidy_metadata$comb_id)
  
  # calc nonzero occurences
  # obj$data$tax_occ <- calc_n_samples(obj, "tax_data", groups = tidy_metadata$PARTICIPANT_ID, 
  #                                    cols = tidy_metadata$comb_id)
  
  # plot heat trees for sample
  # pdf(file.path(graph_dir, "Danish_Diff_Heat_Tree.pdf"))
  # heat_tree(obj,
  #           node_label = taxon_names,
  #           node_size = n_obs,
  #           node_color = obj$data$type_abund$`1_20`,
  #           node_color_axis_label = "Relative Abundance",
  #           node_size_axis_label = "Reads",
  #           layout = "da", # davidson-harel
  #           initial_layout = "re") # fruchterman-reingold
  # dev.off()
  
  for (env_var in env_vars) {
    print(env_var)
    if (env_var == "TRIMESTER") continue # only use bipartite metadata
    
    # compare
    obj$data$diff_table = compare_groups(obj, "type_abund", cols = tidy_metadata$comb_id, groups = tidy_metadata[,env_var])
    print(obj$data$diff_table)
    
    # save text sizes before fdr
    text_size <- obj$n_subtaxa() * 5
    
    # plot heat trees per group
    pdf(file.path(graph_dir, paste0("Danish_Heat_Tree_", env_var, ".pdf")))
    heat_tree(obj,
              node_label = taxon_names,
              node_size = text_size,
              node_color = -median_diff, # we actually want treatment 2 - treatment 1
              node_color_axis_label = "Median difference",
              node_color_interval = c(-10, 10), # symmetric
              node_color_range = c("cyan", "gray", "magenta"),
              layout = "da", # davidson-harel
              initial_layout = "re" # fruchterman-reingold
    ) %>% print()
    dev.off()
    
    # write data
    write.csv(obj$data$diff_table, file.path(graph_dir, paste0("Danish_Taxon_Data_", env_var, ".csv")))
    
    # fdr correct
    obj <- mutate_obs(obj, "diff_table", wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"))
    obj$data$diff_table$log2_median_ratio[obj$data$diff_table$wilcox_p_value > 0.05] <- 0 # filter insignificant
    pdf(file.path(graph_dir, paste0("Danish_Heat_Tree_Significant_", env_var, ".pdf")))
    heat_tree(obj,
              node_label = taxon_names,
              node_size = text_size,
              node_color = -median_diff, # we actually want treatment 2 - treatment 1
              node_color_axis_label = "Median difference",
              node_color_interval = c(-10, 10), # symmetric
              node_color_range = c("cyan", "gray", "magenta"),
              layout = "da", # davidson-harel
              initial_layout = "re" # fruchterman-reingold
    ) %>% print()
    dev.off()
  }
}

diff_heat_trees("Metaphlan")
