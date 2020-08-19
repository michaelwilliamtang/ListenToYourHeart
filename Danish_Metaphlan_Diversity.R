# 8/17/20
# Calculates alpha and beta diversity for metaphlan at genus level

library(tidyverse)
library(vegan)
data_dir <- file.path("Data", "Tidy_Danish")
summarize <- dplyr::summarize

# load
ds1 <- "Metaphlan_Genus"
load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))

# spread into matrix
tidy_df <- tidy_df %>% 
  spread(analyte, val)
tidy_mat <- tidy_df %>%
  select(-comb_id) %>%
  as.matrix()
# rownames(tidy_mat) = tidy_df$comb_id

# calc alpha diversity indices
shannon_vec <- diversity(tidy_mat, index = "shannon")
simpson_vec <- diversity(tidy_mat, index = "simpson")

# calc beta diversity indices
# beta_vec <- betadiver(tidy_mat, "w")

# gather and save
tidy_df <- cbind.data.frame(comb_id = tidy_df$comb_id, shannon = shannon_vec, simpson = simpson_vec)
tidy_df <- tidy_df %>% 
  gather("analyte", "val", -comb_id)
save(tidy_df, tidy_metadata, file = file.path(data_dir, "Tidy_Metaphlan_Diversity.RData"))
