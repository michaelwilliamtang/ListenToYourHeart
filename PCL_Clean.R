# 6/10/20
# Tidies and cleans PCL data, combining batched/unbatched/omni/regular

library(tidyverse)
# library(fs)
data_dir <- "Data"
save_dir <- "Data/Tidy"
summarize <- dplyr::summarize

# get unbatched data (no participant ids, so we must match those from batched data)
pcl_df <- read.table(file.path(data_dir, "pathabundance_april_30.csv"), sep = ",", stringsAsFactors = FALSE)
pcl_df$V1[50] <- paste(pcl_df$V1[50], "b", sep = "_") # dup row
full_names <- pcl_df$V1
full_names[184:length(full_names)] <- full_names[184:length(full_names)] %>% 
  str_replace_all("[\\||:|\\s|\\-|_]", "\\.")
# full_names <- pcl_df$V1 %>% path_sanitize()
rep_names <- full_names
rep_names[184:length(full_names)] <- 184:length(full_names)
rownames(pcl_df) <- rep_names
pcl_df <- pcl_df %>%
  select(-V1) %>%
  as.matrix() %>% t() %>%
  as_tibble()
# pcl_df %>% select(study_id, participant_id) %>% View()
# View(pcl_df[,1:100])
pcl_df$study_id <- str_replace_all(pcl_df$study_id, " ", "")
pathways <- grep("\\.\\.", full_names, perl = T)
pcl_metadata <- pcl_df %>% select(which(!grepl("\\.\\.", full_names, perl = T)))
pcl_pathway_df <- pcl_df %>%
  select(study_id, collection_type, omni_comparison, pathways)
# pcl_pathway_df$participant_id[16] <- paste(pcl_pathway_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# pcl_pathway_df$participant_id[18] <- paste(pcl_pathway_df$participant_id[18], "b", sep = "_")
pcl_pathway_df <- pcl_pathway_df %>%
  gather("pathway", "val", -c(study_id, collection_type, omni_comparison))
# pcl_pathway_df$pathway <- str_replace_all(pcl_pathway_df$pathway, "[//|| |:]", "\\.")
pcl_pathway_df$val <- as.numeric(pcl_pathway_df$val)
pcl_pathway_df$val <- log2(pcl_pathway_df$val + 1) # allows for 0's
pcl_pathway_df$omni_comparison[which(pcl_pathway_df$omni_comparison == "")] <- "no"
pcl_pathway_df$pathway <- as.numeric(pcl_pathway_df$pathway)
pcl_pathway_df <- pcl_pathway_df %>%
  mutate(pathway = full_names[pathway])

# get batched data
pcl_df <- read.table(file.path(data_dir, "hf-abundance-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
pcl_df <- pcl_df %>%
  as.matrix() %>% t() %>%
  as_tibble()
# View(pcl_df[,1:100])
pathways <- grep("\\.\\.", colnames(pcl_df), perl = T)
pcl_pathway_df2 <- pcl_df %>%
  select(study_id, participant_id, collection_type, omni_comparison, pathways)
# pcl_pathway_df2$participant_id[16] <- paste(pcl_pathway_df2$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# pcl_pathway_df2$participant_id[18] <- paste(pcl_pathway_df2$participant_id[18], "b", sep = "_")
pcl_pathway_df2 <- pcl_pathway_df2 %>%
  gather("pathway", "val", -c(study_id, participant_id, omni_comparison, collection_type))
pcl_pathway_df2$val <- as.numeric(pcl_pathway_df2$val)
pcl_pathway_df2 <- pcl_pathway_df2 %>%
  mutate(batch_val = val)

# take only samples that have both batched and unbatched versions
unbatched_ids <- pcl_pathway_df$study_id %>% unique()
batched_ids <- pcl_pathway_df2$study_id %>% unique()
valid_ids <- unbatched_ids[which(unbatched_ids %in% batched_ids)]
pcl_pathway_df <- pcl_pathway_df %>% filter(study_id %in% valid_ids)
pcl_pathway_df2 <- pcl_pathway_df2 %>% filter(study_id %in% valid_ids)

# take only pathways that have both batched and unbatched versions
unbatched_paths <- pcl_pathway_df$pathway %>% unique()
batched_paths <- pcl_pathway_df2$pathway %>% unique()
valid_paths <- unbatched_paths[which(unbatched_paths %in% batched_paths)]
pcl_pathway_df <- pcl_pathway_df %>% filter(pathway %in% valid_paths)
pcl_pathway_df2 <- pcl_pathway_df2 %>% filter(pathway %in% valid_paths)

# check before merge
pcl_pathway_df <- pcl_pathway_df %>% arrange(pathway, study_id, collection_type)
pcl_pathway_df2 <- pcl_pathway_df2 %>% arrange(pathway, study_id, collection_type)
table(pcl_pathway_df$study_id == pcl_pathway_df2$study_id)
table(pcl_pathway_df$pathway == pcl_pathway_df2$pathway)
table(pcl_pathway_df$collection_type == pcl_pathway_df2$collection_type)
# table(pcl_pathway_df$omni_comparison == pcl_pathway_df2$omni_comparison)

# merge
pcl_pathway_df2$omni_comparison <- pcl_pathway_df$omni_comparison # inconsistent, so we replace one set with the other
pcl_pathway_df <- cbind(participant_id = pcl_pathway_df2$participant_id, pcl_pathway_df) # unbatched data has no part ids
comp_batch_pathway <- cbind(pcl_pathway_df, batch_val = pcl_pathway_df2$batch_val)
pcl_pathway_df$batch <- F
pcl_pathway_df2$batch <- T
pcl_pathway_df2$val <- pcl_pathway_df2$batch_val
pcl_pathway_df <- rbind(pcl_pathway_df, pcl_pathway_df2 %>% select(-batch_val))

save(pcl_pathway_df, comp_batch_pathway, pcl_metadata, file = file.path(save_dir, "Tidy_PCL.RData"))