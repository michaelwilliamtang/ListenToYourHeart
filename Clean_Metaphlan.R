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
taxa <- 183:ncol(metaphlan_df)
# met_metadata <- metaphlan_df %>% select(1:182)
met_taxa_df <- metaphlan_df %>%
  select(study_id, collection_type, omni_comparison, taxa)
# unbatched_spread_df <- met_taxa_df[4:ncol(met_taxa_df)]
# unbatched_meta_df <- met_taxa_df[1:3]
# met_taxa_df$participant_id[16] <- paste(met_taxa_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# met_taxa_df$participant_id[18] <- paste(met_taxa_df$participant_id[18], "b", sep = "_")
met_taxa_df <- met_taxa_df %>%
  gather("taxa", "val", -c(study_id, collection_type, omni_comparison))
met_taxa_df$taxa <- str_replace_all(met_taxa_df$taxa, "\\|", "\\.")
met_taxa_df$val <- as.numeric(met_taxa_df$val)
met_taxa_df$val <- log2(met_taxa_df$val + 1) # allows for 0's
met_taxa_df$omni_comparison[which(met_taxa_df$omni_comparison == "")] <- "no"
# unbatched_meta_df$omni_comparison[which(unbatched_meta_df$omni_comparison == "")] <- "no"
# clean spread df too
# colnames(unbatched_spread_df) <- colnames(unbatched_spread_df) %>%
#   str_replace_all("\\|", "\\.")
# unbatched_spread_df <- unbatched_spread_df %>% 
#   sapply(as.numeric)
# unbatched_spread_df <- log2(unbatched_spread_df + 1)


# get batched data
# read metaphlan
metaphlan_df <- read.table(file.path(data_dir, "hf-metaphlan-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
metaphlan_df <- metaphlan_df %>%
  as.matrix() %>% t() %>%
  as_tibble()
# View(metaphlan_df[,1:100])
taxa <- 183:ncol(metaphlan_df)
met_metadata <- metaphlan_df %>% select(1:182)
met_metadata$omni_comparison[which(met_metadata$omni_comparison == "")] <- "no"
met_taxa_df2 <- metaphlan_df %>%
  select(study_id, participant_id, collection_type, omni_comparison, taxa)
# unbatched_part_ids <- met_taxa_df2$participant_id # since unbatched data did not have any, must borrow
# names(unbatched_part_ids) <- met_taxa_df2$study_id
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

# add taxum level info
t_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
met_taxa_df$level <- t_levels[str_count(met_taxa_df$taxa, "__")]

# add part ids to unbatched
# unbatched_meta_df <- unbatched_meta_df %>%
  # mutate(participant_id = unbatched_part_ids[study_id])

# filter lvad samples
lvad <- met_metadata$lvad
names(lvad) <- met_metadata$study_id
met_taxa_df <- met_taxa_df %>% filter(lvad[study_id] == 0)

# save
# tmp <- unbatched_spread_df
# unbatched_spread_df <- cbind(unbatched_meta_df, tmp)
save(met_taxa_df, met_metadata, file = file.path(save_dir, "Tidy_Log_Metaphlan.RData"))

met_taxa_df$level <- t_levels[str_count(met_taxa_df$taxa, "__")]
met_taxa_df$val <- 2 ^ met_taxa_df$val - 1
# shift by max neg
min_val <- min(met_taxa_df$val)
if (min_val < 0) met_taxa_df <- met_taxa_df %>% mutate(val = val - min_val)

# tmp <- 2 ^ tmp - 1
# unbatched_spread_df <- cbind(unbatched_meta_df, tmp)
save(met_taxa_df, met_metadata, file = file.path(save_dir, "Tidy_Metaphlan.RData"))

met_levs <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
met_taxa_df <- met_taxa_df %>% group_by(batch, level, study_id) %>%
  mutate(total = sum(val)) %>%
  mutate(max = max(total)) %>%
  mutate(val = val / max) %>%
  ungroup() %>%
  select(-total, -max)
# met_taxa_df$val <- met_taxa_df$val / 100
# tmp <- tmp / 100
# unbatched_spread_df <- cbind(unbatched_meta_df, tmp)
save(met_taxa_df, met_metadata, file = file.path(save_dir, "Tidy_Scaled_Metaphlan.RData"))