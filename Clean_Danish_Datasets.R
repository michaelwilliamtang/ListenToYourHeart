# 6/25/20
# Tidies and cleans metaphlan data, combining batched/unbatched/omni/regular

library(tidyverse)
data_dir <- file.path("Data", "Danish")
save_dir <- file.path("Data", "Tidy_Danish")
if (!dir.exists(save_dir)) dir.create(save_dir)
summarize <- dplyr::summarize

### cleaning metaphlan genus
# read
ds <- "Metaphlan_Genus"
load(file.path(data_dir, paste("Full_Danish", ds, "df.RData", sep = "_")))
met_taxa_df <- pcl_df2
met_metadata <- pcl_metadata
met_metadata <- met_metadata %>% t() %>% as.data.frame()
met_metadata[] <- lapply(met_metadata, as.character)
table(rownames(met_metadata) == rownames(met_taxa_df)) # check match

# match, build comb id
met_taxa_df <- cbind(comb_id = paste(met_metadata$PARTICIPANT_ID, met_metadata$GESTATIONAL_WEEK, sep = "_"), met_taxa_df,
                     stringsAsFactors = F)
table(met_taxa_df$comb_id %>% duplicated()) # check dup
met_taxa_df <- met_taxa_df %>%
  gather("taxa", "val", -comb_id)

# clean name to genus level
met_reg <- ".*g__"
met_taxa_df <- met_taxa_df %>% mutate(taxa = str_replace(taxa, met_reg, ""))
save(met_taxa_df, met_metadata, file = file.path(save_dir, paste0("Tidy_", ds, ".RData")))

### cleaning uBiome
# read
ds <- "uBiome"
load(file.path(data_dir, paste("Full_Danish", ds, "df.RData", sep = "_")))
ub_taxa_df <- pcl_df2
ub_metadata <- pcl_metadata
ub_metadata <- ub_metadata %>% t() %>% as.data.frame()
ub_metadata[] <- lapply(ub_metadata, as.character)
table(rownames(ub_metadata) == rownames(ub_taxa_df)) # check match

# match, build comb id
# ub_taxa_df <- cbind(comb_id = paste(ub_metadata$PARTICIPANT_ID, ub_metadata$GESTATIONAL_WEEK, sep = "_"), ub_taxa_df,
#                      stringsAsFactors = F)
ub_taxa_df <- cbind(comb_id = rownames(ub_taxa_df), ub_taxa_df,
                   stringsAsFactors = F)
table(ub_taxa_df$comb_id %>% duplicated()) # check dup
ub_taxa_df <- ub_taxa_df %>%
  gather("taxa", "val", -comb_id)
ub_taxa_df$val[which(is.na(ub_taxa_df$val))] = 0

# ub_reg <- ".*g__"
# ub_taxa_df <- ub_taxa_df %>% mutate(taxa = str_replace(taxa, ub_reg, ""))
save(ub_taxa_df, ub_metadata, file = file.path(save_dir, paste0("Tidy_", ds, ".RData")))

### cleaning PCL
# read
ds <- "pathabundance"
ds2 <- "PCL"
load(file.path(data_dir, paste("Full_Danish", ds, "df.RData", sep = "_")))
pcl_pathway_df <- pcl_df2
pcl_metadata <- pcl_metadata
pcl_metadata <- pcl_metadata %>% t() %>% as.data.frame()
pcl_metadata[] <- lapply(pcl_metadata, as.character)
table(rownames(pcl_metadata) == rownames(pcl_pathway_df)) # check match

# match, build comb id
pcl_pathway_df <- cbind(comb_id = paste(pcl_metadata$PARTICIPANT_ID, pcl_metadata$GESTATIONAL_WEEK, sep = "_"), pcl_pathway_df,
                     stringsAsFactors = F)
table(pcl_pathway_df$comb_id %>% duplicated()) # check dup
pcl_pathway_df <- pcl_pathway_df %>%
  gather("pathway", "val", -comb_id)

# ub_reg <- ".*g__"
# ub_taxa_df <- ub_taxa_df %>% mutate(taxa = str_replace(taxa, ub_reg, ""))
save(ub_taxa_df, ub_metadata, file = file.path(save_dir, paste0("Tidy_", ds2, ".RData")))

### cleaning full metaphlan
# read
ds <- "Metaphlan"
metaphlan_df <- read.table(file.path(data_dir, "metaphlan.pcl"), sep = "\t", row.names = 1, stringsAsFactors = F)
met_taxa_df <- metaphlan_df %>% t() %>% as.data.frame(stringsAsFactors = F)
colnames(met_taxa_df) <- colnames(met_taxa_df) %>% str_replace_all(" ", "_")
colnames(met_taxa_df) <- colnames(met_taxa_df) %>% str_replace_all("\\(|\\)", "\\.")

# match, build comb id
met_taxa_df <- met_taxa_df %>% mutate(comb_id = paste(PARTICIPANT_ID, GESTATIONAL_WEEK, sep = "_"))
table(met_taxa_df$comb_id %>% duplicated()) # check dup
met_taxa_df <- met_taxa_df %>%
  select(-c(SampleID, ORIGINAL_SAMPLE_ID, GESTATIONAL_WEEK, TRIMESTER, MOTHER_AGE_AT_CONCEPTION_.DAYS., 
            BABY_BIRTHWEIGHT_.G., BABY_LENGTH_.CM., PARTICIPANT_ID)) %>%
  gather("taxa", "val", -comb_id)

save(met_taxa_df, file = file.path(save_dir, paste0("Tidy_", ds, ".RData")))

### cleaning meta pathcoverage
# read
ds <- "Meta_Pathcoverage"
metaphlan_df <- read.table(file.path(data_dir, "meta-pathcoverage.pcl"), sep = "\t", row.names = 1, stringsAsFactors = F)
met_taxa_df <- metaphlan_df %>% t() %>% as.data.frame(stringsAsFactors = F)
colnames(met_taxa_df) <- colnames(met_taxa_df) %>% str_replace_all(" ", "_")
colnames(met_taxa_df) <- colnames(met_taxa_df) %>% str_replace_all("\\(|\\)", "\\.")
full_names <- colnames(met_taxa_df)
rep_names <- full_names
rep_names[9:length(rep_names)] = 9:length(rep_names)
colnames(met_taxa_df) <- rep_names

# match, build comb id
met_taxa_df <- met_taxa_df %>% mutate(comb_id = paste(PARTICIPANT_ID, GESTATIONAL_WEEK, sep = "_"))
table(met_taxa_df$comb_id %>% duplicated()) # check dup
met_taxa_df <- met_taxa_df %>%
  select(-c(PATHWAY, ORIGINAL_SAMPLE_ID, GESTATIONAL_WEEK, TRIMESTER, MOTHER_AGE_AT_CONCEPTION_.DAYS., 
            BABY_BIRTHWEIGHT_.G., BABY_LENGTH_.CM., PARTICIPANT_ID)) %>%
  gather("taxa", "val", -comb_id)
met_taxa_df$taxa <- met_taxa_df$taxa %>% as.numeric()
met_taxa_df$taxa <- full_names[met_taxa_df$taxa]

save(met_taxa_df, file = file.path(save_dir, paste0("Tidy_", ds, ".RData")))