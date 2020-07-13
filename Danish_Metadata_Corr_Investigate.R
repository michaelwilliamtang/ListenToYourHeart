# 7/13/20
# Filter and plot analyte-metadatum pairs from Danish data

library(tidyverse)
data_dir <- "Data/Danish/Danish_From_GCloud"
graph_dir <- "Graphs/Danish"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize

fdr_threshold <- 0.05

danish_corr_clean <- function() {
  # file_name <- "uBiomeStatsDanishPart.txt"
  file_name <- "pclStatsDanishPart.txt"
  file_name <- "metaphlan/pclStatsDanishPart.txt"
  data_df <- read.table(file.path(data_dir, file_name), sep = "\t", row.names = 1, stringsAsFactors = F, 
                         quote = "")
 
}

# gestational weeks, mother
clean_names <- c("Gestational Age (weeks)", "Mother's Age at Birth (years)", 
                 "Baby Birth Weight (g)", "Baby Birth Length (cm)")
analytes <- c("X.GW.", "X.mage.", "X.BB.", "X.BL.")
analytes2 <- c("X.GW_fdr.", "X.mage_fdr.", "X.BB_fdr.", "X.BL_fdr.")