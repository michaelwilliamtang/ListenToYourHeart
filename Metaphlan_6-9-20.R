library(tidyverse)
data_dir <- "Data"
graph_dir <- "Graphs/Metaphlan"
summarize <- dplyr::summarize
desc <- "Log_Unbatch"

# read metaphlan
metaphlan_df <- read.table(file.path(data_dir, "metaphlan2_april_30.csv"), sep = ",", stringsAsFactors = FALSE)
metaphlan_df$V1[50] <- paste(metaphlan_df$V1[50], "b", sep = "_") # dup row
rownames(metaphlan_df) <- metaphlan_df$V1
metaphlan_df <- metaphlan_df %>%
  select(-V1) %>%
  as.matrix() %>% t() %>%
  as_tibble()
metaphlan_df$participant_id <- metaphlan_df$study_id %>% str_sub(1,6)
# metaphlan_df %>% select(study_id, participant_id) %>% View()
# View(metaphlan_df[,1:100])

# get top 5 phyla
phyla <- grep(".*p__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_phyla_df <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, phyla)
colnames(met_phyla_df) <- met_phyla_df %>% colnames() %>%
  str_replace(".*p__", "")
met_phyla_df$participant_id[16] <- paste(met_phyla_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
met_phyla_df$participant_id[18] <- paste(met_phyla_df$participant_id[18], "b", sep = "_")
met_phyla_df <- met_phyla_df %>%
  gather("taxa", "val", -c(participant_id, collection_type))
met_phyla_df$val <- as.numeric(met_phyla_df$val)
met_phyla_df$val <- log2(met_phyla_df$val + 1) # allows for 0's
phy <- met_phyla_df %>%
  group_by(taxa) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_5_phy <- phy$taxa[1:5]

# plot most diff phyla, omni vs reg
comp_phyla <- met_phyla_df %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  spread(collection_type, mean_val) %>%
  mutate(collection_diff = omni - regular,
         collection_diff_abs = abs(collection_diff)) %>%
  arrange(desc(collection_diff_abs))
diff_20_phy <- comp_phyla$taxa[1:20]
comp_phyla <- comp_phyla %>% filter(taxa %in% diff_20_phy)
pdf(file.path(graph_dir, paste(desc, "Phyla_By_Collection_Diff.pdf", sep = "_")), width = 6, height = 4)
comp_phyla %>% ggplot(aes(x = reorder(taxa, -collection_diff), y = collection_diff)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Phylum") +
  ylab("Relative Abundance Diff (omni - reg)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# plot top phyla
met_phyla_df2 <- met_phyla_df %>%
  filter(taxa %in% top_5_phy) %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(met_phyla_df)
pdf(file.path(graph_dir, paste(desc, "Top_5_Phyla_Side.pdf", sep = "_")), width = 6, height = 4)
met_phyla_df2 %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Phylum") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pdf(file.path(graph_dir, paste(desc, "Top_5_Phyla_Stack.pdf", sep = "_")), width = 6, height = 4)
met_phyla_df2 %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Phylum") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
met_phyla_df2 <- met_phyla_df %>%
  filter(taxa %in% top_5_phy) %>%
  group_by(collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
pdf(file.path(graph_dir, paste(desc, "Phyla_Aggregated.pdf", sep = "_")), width = 6, height = 4)
met_phyla_df2 %>% ggplot(aes(x = reorder(collection_type, -mean_val), y = mean_val)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Collection Type") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# get top 20 genera
genera <- grep(".*g__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_genera_df <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, genera)
colnames(met_genera_df) <- met_genera_df %>% colnames() %>%
  str_replace(".*g__", "")
met_genera_df$participant_id[16] <- paste(met_genera_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
met_genera_df$participant_id[18] <- paste(met_genera_df$participant_id[18], "b", sep = "_")
met_genera_df <- met_genera_df %>%
  gather("taxa", "val", -c(participant_id, collection_type))
met_genera_df$val <- as.numeric(met_genera_df$val)
met_genera_df$val <- log2(met_genera_df$val + 1) # allows for 0's
gen <- met_genera_df %>%
  group_by(taxa) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_5_gen <- gen$taxa[1:20]

# plot most diff genera, omni vs reg
comp_genera <- met_genera_df %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  spread(collection_type, mean_val) %>%
  mutate(collection_diff = omni - regular,
         collection_diff_abs = abs(collection_diff)) %>%
  arrange(desc(collection_diff_abs))
diff_20_gen <- comp_genera$taxa[1:20]
comp_genera <- comp_genera %>% filter(taxa %in% diff_20_gen)
pdf(file.path(graph_dir, paste(desc, "Genera_By_Collection_Diff.pdf", sep = "_")), width = 6, height = 4)
comp_genera %>% ggplot(aes(x = reorder(taxa, -collection_diff), y = collection_diff)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Genus") +
  ylab("Relative Abundance Diff (omni - reg)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# plot top genera
met_genera_df2 <- met_genera_df %>%
  filter(taxa %in% top_5_gen) %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(met_genera_df)
pdf(file.path(graph_dir, paste(desc, "Top_20_Genera_Side.pdf", sep = "_")), width = 6, height = 4)
met_genera_df2 %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Genus") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pdf(file.path(graph_dir, paste(desc, "Top_20_Genera_Stack.pdf", sep = "_")), width = 6, height = 4)
met_genera_df2 %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Genus") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
met_genera_df2 <- met_genera_df %>%
  filter(taxa %in% top_5_gen) %>%
  group_by(collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
pdf(file.path(graph_dir, paste(desc, "Genera_Aggregated.pdf", sep = "_")), width = 6, height = 4)
met_genera_df2 %>% ggplot(aes(x = reorder(collection_type, -mean_val), y = mean_val)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Collection Type") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

