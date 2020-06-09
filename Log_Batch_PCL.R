# 6/8/20
# Compares omni and normal for log, batch-corrected PCL

library(tidyverse)
data_dir <- "Data"
graph_dir <- "Graphs/PCL"
summarize <- dplyr::summarize

# read pcl
pcl_df <- read.table(file.path(data_dir, "hf-abundance-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
pcl_df <- pcl_df %>% 
  as.matrix() %>% t() %>% 
  as_tibble()
# View(pcl_df[,1:50])

# get top 40 pathways
pathways <- grep("\\.\\.", colnames(pcl_df), perl = T)
pcl_pathways_df <- pcl_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, pathways)
pcl_pathways_df$participant_id[16] <- paste(pcl_pathways_df$participant_id[16], "b", sep = "") # anomaly participant (dup)
pcl_pathways_df$participant_id[18] <- paste(pcl_pathways_df$participant_id[18], "b", sep = "")
pcl_pathways_df <- pcl_pathways_df %>%
  gather("pathway", "val", -c(participant_id, collection_type))
pcl_pathways_df$val <- as.numeric(pcl_pathways_df$val)
path <- pcl_pathways_df %>%
  group_by(pathway) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_40_path <- path$pathway[1:40]

# plot most diff pathways, omni vs reg
comp_pathways <- pcl_pathways_df %>%
  group_by(pathway, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  spread(collection_type, mean_val) %>%
  mutate(collection_diff = omni - regular,
         collection_diff_abs = abs(collection_diff)) %>%
  arrange(desc(collection_diff_abs))
diff_20_path <- comp_pathways$pathway[1:20]
comp_pathways <- comp_pathways %>% filter(pathway %in% diff_20_path)
pdf(file.path(graph_dir, "Pathways_By_Collection_Diff.pdf"), width = 18, height = 12)
comp_pathways %>% ggplot(aes(x = reorder(pathway, -collection_diff), y = collection_diff)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance Diff (omni - reg)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# plot top pathways
pcl_pathways_df2 <- pcl_pathways_df %>%
  filter(pathway %in% top_40_path) %>%
  group_by(pathway, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(pcl_pathways_df)
pdf(file.path(graph_dir, "Top_40_Pathways_Side.pdf"), width = 18, height = 12)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(pathway, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pdf(file.path(graph_dir, "Top_40_Pathways_Stack.pdf"), width = 18, height = 12)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(pathway, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pcl_pathways_df2 <- pcl_pathways_df %>%
  filter(pathway %in% top_40_path) %>%
  group_by(collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
pdf(file.path(graph_dir, "Pathways_Aggregated.pdf"), width = 6, height = 4)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(collection_type, -mean_val), y = mean_val)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Collection Type") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()