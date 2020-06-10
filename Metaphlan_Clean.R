# 6/10/20
# Tidies and cleans metaphlan data, combining batched/unbatched/omni/regular

library(tidyverse)
data_dir <- "Data"
save_dir <- "Data/Tidy"
summarize <- dplyr::summarize

# get unbatched data (no participant ids, so we must match those from batched data)
metaphlan_df <- read.table(file.path(data_dir, "metaphlan2_april_30.csv"), sep = ",", stringsAsFactors = FALSE)
metaphlan_df$V1[50] <- paste(metaphlan_df$V1[50], "b", sep = "_") # dup row
rownames(metaphlan_df) <- metaphlan_df$V1
metaphlan_df <- metaphlan_df %>%
  select(-V1) %>%
  as.matrix() %>% t() %>%
  as_tibble()
# metaphlan_df %>% select(study_id, participant_id) %>% View()
# View(metaphlan_df[,1:100])
metaphlan_df$study_id <- str_replace_all(metaphlan_df$study_id, " ", "")
taxa <- grep("__", colnames(metaphlan_df), perl = T)
met_metadata <- metaphlan_df %>% select(which(!grepl(".*__[A-Za-z]+$", colnames(metaphlan_df), perl = T)))
met_taxa_df <- metaphlan_df %>%
  select(study_id, collection_type, omni_comparison, taxa)
# met_taxa_df$participant_id[16] <- paste(met_taxa_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# met_taxa_df$participant_id[18] <- paste(met_taxa_df$participant_id[18], "b", sep = "_")
met_taxa_df <- met_taxa_df %>%
  gather("taxa", "val", -c(study_id, collection_type, omni_comparison))
met_taxa_df$taxa <- str_replace_all(met_taxa_df$taxa, "\\|", "\\.")
met_taxa_df$val <- as.numeric(met_taxa_df$val)
met_taxa_df$val <- log2(met_taxa_df$val + 1) # allows for 0's
met_taxa_df$omni_comparison[which(met_taxa_df$omni_comparison == "")] <- "no"

# get batched data
# read metaphlan
metaphlan_df <- read.table(file.path(data_dir, "hf-metaphlan-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
metaphlan_df <- metaphlan_df %>%
  as.matrix() %>% t() %>%
  as_tibble()
# View(metaphlan_df[,1:100])
taxa <- grep("__", colnames(metaphlan_df), perl = T)
met_taxa_df2 <- metaphlan_df %>%
  select(study_id, participant_id, collection_type, omni_comparison, taxa)
# met_taxa_df2$participant_id[16] <- paste(met_taxa_df2$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# met_taxa_df2$participant_id[18] <- paste(met_taxa_df2$participant_id[18], "b", sep = "_")
met_taxa_df2 <- met_taxa_df2 %>%
  gather("taxa", "val", -c(study_id, participant_id, omni_comparison, collection_type))
met_taxa_df2$val <- as.numeric(met_taxa_df2$val)
met_taxa_df2 <- met_taxa_df2 %>%
  mutate(batch_val = val)

# take only samples that have both batched and unbatched versions
unbatched_ids <- met_taxa_df$study_id %>% unique()
batched_ids <- met_taxa_df2$study_id %>% unique()
valid_ids <- unbatched_ids[which(unbatched_ids %in% batched_ids)]
met_taxa_df <- met_taxa_df %>% filter(study_id %in% valid_ids)
met_taxa_df2 <- met_taxa_df2 %>% filter(study_id %in% valid_ids)

# take only taxa that have both batched and unbatched versions
unbatched_taxa <- met_taxa_df$taxa %>% unique()
batched_taxa <- met_taxa_df2$taxa %>% unique()
valid_taxa <- unbatched_taxa[which(unbatched_taxa %in% batched_taxa)]
met_taxa_df <- met_taxa_df %>% filter(taxa %in% valid_taxa)
met_taxa_df2 <- met_taxa_df2 %>% filter(taxa %in% valid_taxa)

# check before merge
table(met_taxa_df$study_id == met_taxa_df2$study_id)
table(met_taxa_df$taxa == met_taxa_df2$taxa)
table(met_taxa_df$collection_type == met_taxa_df2$collection_type)
table(met_taxa_df$omni_comparison == met_taxa_df2$omni_comparison)

# merge
met_taxa_df <- cbind(participant_id = met_taxa_df2$participant_id, met_taxa_df) # unbatched data has no part ids
comp_batch_taxa <- cbind(met_taxa_df, batch_val = met_taxa_df2$batch_val)
met_taxa_df$batch <- F
met_taxa_df2$batch <- T
met_taxa_df2$val <- met_taxa_df2$batch_val
met_taxa_df <- rbind(met_taxa_df, met_taxa_df2 %>% select(-batch_val))

save(met_taxa_df, comp_batch_taxa, met_metadata, file = file.path(save_dir, "Tidy_Metaphlan.RData"))