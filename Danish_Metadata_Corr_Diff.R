# 8/4/20
# Filter and plot analyte-metadatum pairs from Danish data, summarized as a single heatmap

library(tidyverse)
corr_dir <- file.path("Data", "Danish", "Danish_From_GCloud")
data_dir <- file.path("Data", "Tidy_Danish")
graph_dir <- file.path("Graphs", "Danish")
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize

library(gplots)
palette <- cm.colors

fdr_threshold <- 0.05

danish_corr_diff <- function(comp_N) {
  # ref
  clean_names <- c("Gestational Age (weeks)", "Mother's Age at Birth (days)",
                   "Baby Birth Weight (g)", "Baby Birth Length (cm)")
  analytes2 <- c("X.GW_fdr.", "X.mage_fdr.", "X.BB_fdr.", "X.BL_fdr.") # fdr column name in corr df

  ## metaphlan_genus
  # read and clean
  file_name <- file.path("metaphlan", "pclStatsDanishPart.txt")
  ds1 <- "Metaphlan_Genus"
  corr_df <- read.table(file.path(corr_dir, file_name), sep = "\t", row.names = 1, stringsAsFactors = F,
                         quote = "")
  met_reg <- ".*g__" # metaphlan corr analytes must be regex renamed in the same way as orig data
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace(X.analyte., met_reg, ""))
  corr_df <- corr_df %>% mutate(X.analyte. = str_replace_all(X.analyte., "\"", "")) # remove unwanted quotes
  analytes <- c("GESTATIONAL_WEEK", "MOTHER_AGE_AT_CONCEPTION_.DAYS.", "BABY_BIRTHWEIGHT_.G.", "BABY_LENGTH_.CM.") # column name in metadata
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))

  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr_Diff", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)

  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    dcd_metadatum(comp_N, corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir2)
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

  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr_Diff", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)

  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    dcd_metadatum(comp_N, corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir2)
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
  
  graph_dir2 <- file.path(graph_dir, paste(ds1, "Corr_Diff", sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # plot subsets
  for (i in 1:length(analytes)) {
    meta1 <- analytes[i]
    dcd_metadatum(comp_N, corr_df, tidy_df, tidy_metadata, clean_names[i], meta1, analytes2[i], graph_dir2)
  }
}

dcd_metadatum <- function(comp_N, corr_df, data_df, metadata_df, cname, meta1, meta2, graph_dir3) {
  
  print(meta1)
  
  # convert metadata (default factor for some reason)
  metadata_df[,meta1] <- as.numeric(metadata_df[,meta1])
  
  # get significant
  corr_df <- corr_df %>% filter(!! sym(meta2) < fdr_threshold)
  selected <- corr_df$X.analyte.
  
  # setup
  summ_df <- tibble()
  
  for (sel in selected) {
    tryCatch({
      # get metadatum-analyte pair
      sel_data <- data_df %>% filter(analyte == sel) %>% select(comb_id, val)
      sel_meta <- metadata_df %>% select(comb_id, !! sym(meta1), PARTICIPANT_ID)
      
      # filter missing metadata, merge
      valid_ids <- sel_meta$comb_id[!is.na(sel_meta[,meta1])]
      sel_data <- sel_data %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      sel_meta <- sel_meta %>% filter(comb_id %in% valid_ids) %>% arrange(comb_id)
      # print(table(sel_data$comb_id == sel_meta$comb_id))
      comb_sel <- cbind(sel_data %>% select(val, comb_id), meta = sel_meta[, meta1], id = sel_meta$PARTICIPANT_ID)
      
      # test sum == 0
      comb_sum <- comb_sel$meta %>% sum()
      if (comb_sum == 0) {
        message(paste0("All zeroes with ", sel, " for ", cname))
        next
      }
      
      # calculate low, high, diff, split based on metadata median
      med_meta <- comb_sel$meta %>% median()
      comb_sel <- comb_sel %>% mutate(half = ifelse(meta <= med_meta, "low", "high"))
      comb_sel <- comb_sel %>% group_by(half) %>%
        summarize(mean_val = mean(val)) %>%
        spread(half, mean_val) %>%
        mutate(diff = high - low,
               diff_abs = abs(diff))
      comb_sel$analyte = sel
      summ_df <- rbind(summ_df, comb_sel)
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
  
  # get top N diffs
  summ_df <- summ_df %>% arrange(desc(diff_abs))
  write.table(summ_df, row.names = F, file = file.path(graph_dir3, paste(meta1, "Diff_Table.tsv", sep = "_")), sep = "\t", quote = FALSE)
  comp_N <- min(comp_N, nrow(summ_df)) # set upper bound
  diff_N_anl <- summ_df$analyte[1:comp_N]
  comp_df <- summ_df %>% filter(analyte %in% diff_N_anl)
  gg <- comp_df %>% ggplot(aes(x = reorder(analyte, -diff), y = diff)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab(meta1) +
    ylab("Diff Between Halves, Split By Metadatum\n(top - bottom)") +
    theme(axis.text.x = element_text(size = 6, angle = 90, hjust = 1))
  plot(gg)
  ggsave(file.path(graph_dir3, paste(meta1, "Diff_Bar.pdf", sep = "_")), plot = last_plot())
  
  # plot together in heatmap
  summ_mat <- as.matrix(summ_df %>% select(-analyte))
  rownames(summ_mat) <- summ_df$analyte
  pdf(file.path(graph_dir3, paste(meta1, "Diff_Heatmap.pdf", sep = "_")))
  heatmap.2(summ_mat, Rowv = F, Colv = F, srtCol = 45, cexCol = 1, cexRow = 0.2,
            col = palette, dendrogram = "none", trace = "none", margins = c(8,15))
  dev.off()
}

# danish_corr_diff(80)
