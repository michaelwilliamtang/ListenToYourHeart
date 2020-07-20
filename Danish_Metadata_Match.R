# 7/15/20
# Matches participant ids to uBiome ids through metadata

library(tidyverse)
save_dir <- file.path("Data", "Tidy_Danish")
if (!dir.exists(save_dir)) dir.create(save_dir)
summarize <- dplyr::summarize

# load part ids
file_name <- "Clinical_info_hipaa_compliant.csv"
data_dir <- file.path("Data", "Danish")
part_df <- read.csv(file.path(data_dir, file_name), quote = "")
part_df <- part_df %>% mutate(matcher = paste(GA_Days, Mom_age_birth, Birth_weight_g, sep = "_"))
pid = part_df$Participant.number
names(pid) = part_df$matcher

# load metadata
ds1 <- "uBiome"
data_dir <- save_dir
load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
tidy_metadata <- tidy_metadata %>% mutate(matcher = paste(GA_Days, Mom_age_birth, Birth_weight_g, sep = "_"))

# attach part ids to metadata
tidy_metadata <- tidy_metadata %>% mutate(PARTICIPANT_ID = pid[matcher]) %>%
  select(-matcher)
tidy_metadata$PARTICIPANT_ID <- as.factor(tidy_metadata$PARTICIPANT_ID)

# save
save(tidy_df, tidy_metadata, file = file.path(save_dir, paste0("Tidy_", ds1, ".RData")))