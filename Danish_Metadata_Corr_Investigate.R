# 7/13/20
# Filter and plot analyte-metadatum pairs from Danish data

library(tidyverse)
corr_dir <- "Data/Danish/Danish_From_GCloud"
data_dir <- "Data/Tidy_Danish"
graph_dir <- "Graphs/Danish"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize

fdr_threshold <- 0.05

danish_corr_inv <- function() {
  # ref
  clean_names <- c("Gestational Age (weeks)", "Mother's Age at Birth (days)", 
                   "Baby Birth Weight (g)", "Baby Birth Length (cm)")
  analytes2 <- c("X.GW_fdr.", "X.mage_fdr.", "X.BB_fdr.", "X.BL_fdr.") # fdr column name in corr df
  
  ## metaphlan_genus
  # read and clean
  file_name <- "metaphlan/pclStatsDanishPart.txt"
  ds1 <- "Metaphlan_Genus"
  corr_df <- read.table(file.path(corr_dir, file_name), sep = "\t", row.names = 1, stringsAsFactors = F, 
                         quote = "")
  met_reg <- ".*g__" # metaphlan corr analytes must be regex renamed in the same way as orig data
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace(X.analyte., met_reg, ""))
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace_all(X.analyte., "\"", "")) # remove unwanted quotes
  analytes <- c("GESTATIONAL_WEEK", "MOTHER_AGE_AT_CONCEPTION_.DAYS.", "BABY_BIRTHWEIGHT_.G.", "BABY_LENGTH_.CM.") # column name in metadata
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    graph_dir3 <- file.path(graph_dir2, meta1)
    if (!dir.exists(graph_dir3)) dir.create(graph_dir3)
    dci_metadatum(corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir3)
  }
  
  
  ## pcl
  # read and clean
  file_name <- "pclStatsDanishPart.txt"
  ds1 <- "PCL"
  corr_df <- read.table(file.path(corr_dir, file_name), sep = "\t", row.names = 1, stringsAsFactors = F, 
                        quote = "")
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace_all(X.analyte., "\"", "")) # remove unwanted quotes
  analytes <- c("GESTATIONAL_WEEK", "MOTHER_AGE_AT_CONCEPTION_.DAYS.", "BABY_BIRTHWEIGHT_.G.", "BABY_LENGTH_.CM.") # column name in metadata
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    graph_dir3 <- file.path(graph_dir2, meta1)
    if (!dir.exists(graph_dir3)) dir.create(graph_dir3)
    dci_metadatum(corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir3)
  }
  
  # diff ref for uBiome
  clean_names <- c("Gestational Age (days)", "Pre-Pregnancy BMI", 
                   "Baby Birth Weight (g)")
  analytes2 <- c("X.GA_Days_fdr.", "X.Pre_Preg_BMI_fdr.", "X.Birth_weight_g_fdr.") # fdr column name in corr df
  
  ## uBiome
  # read and clean
  file_name <- "uBiomeStatsDanishPart.txt"
  ds1 <- "uBiome"
  corr_df <- read.table(file.path(corr_dir, file_name), sep = "\t", row.names = 1, stringsAsFactors = F, 
                        quote = "")
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace_all(X.names., "\"", "")) # remove unwanted quotes
  analytes <- c("GA_Days", "Pre_pregy_weight", "Birth_weight_g") # column name in metadata
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  
  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    graph_dir3 <- file.path(graph_dir2, meta1)
    if (!dir.exists(graph_dir3)) dir.create(graph_dir3)
    dci_metadatum(corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir3)
  }
}

dci_metadatum <- function(corr_df, data_df, metadata_df, cname, meta1, meta2, graph_dir3) {
  
  # convert metadata (default factor for some reason)
  metadata_df[,meta1] <- as.numeric(metadata_df[,meta1])
  
  # get significant
  corr_df <- corr_df %>% filter(!! sym(meta2) < fdr_threshold)
  selected <- corr_df$X.analyte.
  
  for (sel in selected) {
    tryCatch({
      # get metadatum-analyte pair
      sel_data <- data_df %>% filter(analyte == sel) %>% select(comb_id, val)
      sel_meta <- metadata_df %>% select(comb_id, !! sym(meta1))
      
      # filter missing metadata, merge
      valid_ids <- sel_meta$comb_id[!is.na(sel_meta[,meta1])]
      sel_data <- sel_data %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      sel_meta <- sel_meta %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      print(table(sel_data$comb_id == sel_meta$comb_id))
      comb_sel <- cbind(sel_data %>% select(val, comb_id), meta = sel_meta[,meta1])
      
      # plot
      desc <- paste(cname, "vs", sel, sep = " ")
      print(desc)
      gg <- ggplot(comb_sel, aes(x = meta, y = val)) +
        geom_point(size = 2) +
        # geom_text(label = comb_sel$comb_id, 
        #           position = position_dodge(width = 1)) +
        xlab(cname) +
        ylab(sel) +
        labs(title = desc)
      
      plot(gg)
      ggsave(file.path(graph_dir3, paste(sel, meta1, "Corr.pdf", sep = "_")), plot = last_plot())
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

danish_corr_inv()


