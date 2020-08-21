# 8/20/20
# Produces clean, ungathered Danish datasets

library(tidyverse)
data_dir <- file.path("Data", "Tidy_Danish")
save_dir <- file.path("Data", "Tidy_Ungathered_Danish")
if (!dir.exists(data_dir2)) dir.create(save_dir)
summarize <- dplyr::summarize

ungather_dataset <- function(ds) {
  load(file.path(data_dir, paste0("Tidy_", ds, ".RData")))
  tidy_df <- tidy_df %>% 
    spread(analyte, val)
  save(tidy_df, tidy_metadata, file = file.path(save_dir, paste0("Tidy_Ungathered_", ds, ".RData")))
}

datasets <- c("Metaphlan", "uBiome")
for (ds in datasets) ungather_dataset(ds)
