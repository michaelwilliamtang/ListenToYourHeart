library(tidyverse)
data_dir <- "Data"
graph_dir <- "Graphs"
summarize <- dplyr::summarize

# read metaphlan
metaphlan_df <- read.table(file.path(data_dir, "hf-metaphlan-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
metaphlan_df <- metaphlan_df %>%
  as.matrix() %>% t() %>%
  as_tibble()
# View(metaphlan_df[,1:100])

# get top 5 phyla
phyla <- grep(".*p__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_phyla_df <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, phyla)
colnames(met_phyla_df) <- met_phyla_df %>% colnames() %>%
  str_replace(".*p__", "")
met_phyla_df$participant_id[16] <- paste(met_phyla_df$participant_id[16], "b", sep = "") # anomaly participant (dup)
met_phyla_df$participant_id[18] <- paste(met_phyla_df$participant_id[18], "b", sep = "")
met_phyla_df <- met_phyla_df %>%
  gather("taxa", "val", -c(participant_id, collection_type))
met_phyla_df$val <- as.numeric(met_phyla_df$val)
phy <- met_phyla_df %>%
  group_by(taxa) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_5_phy <- phy$taxa[1:5]
# plot phyla
met_phyla_df <- met_phyla_df %>%
  filter(taxa %in% top_5_phy) %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(met_phyla_df)
pdf(file.path(graph_dir, "Top_5_Phyla.pdf"), width = 6, height = 4)
met_phyla_df %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Phylum") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust = 1))
dev.off()

# get top 20 genera
genera <- grep(".*g__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_genera_df <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, genera)
colnames(met_genera_df) <- met_genera_df %>% colnames() %>%
  str_replace(".*g__", "")
met_genera_df$participant_id[16] <- paste(met_genera_df$participant_id[16], "b", sep = "") # anomaly participant (dup)
met_genera_df$participant_id[18] <- paste(met_genera_df$participant_id[18], "b", sep = "")
met_genera_df <- met_genera_df %>%
  gather("taxa", "val", -c(participant_id, collection_type))
met_genera_df$val <- as.numeric(met_genera_df$val)
gen <- met_genera_df %>%
  group_by(taxa) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_5_gen <- gen$taxa[1:20]
# plot genera
met_genera_df <- met_genera_df %>%
  filter(taxa %in% top_5_gen) %>%
  group_by(taxa, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(met_genera_df)
pdf(file.path(graph_dir, "Top_20_Genera.pdf"), width = 6, height = 4)
met_genera_df %>% ggplot(aes(x = reorder(taxa, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Genus") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size=10, angle=45, hjust = 1))
dev.off()

