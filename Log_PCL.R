# 6/9/20
# Compares omni and normal for log, un-batch-corrected PCL

library(tidyverse)
data_dir <- "Data"
graph_dir <- "Graphs/PCL"
summarize <- dplyr::summarize
desc <- "Log_Unbatch"

# read pcl
pcl_df <- read.table(file.path(data_dir, "pathabundance_april_30.csv"), sep = ",", stringsAsFactors = FALSE)
pcl_df$V1[50] <- paste(pcl_df$V1[50], "b", sep = "_") # dup row
full_names <- pcl_df$V1
names(full_names) = full_names
full_names[184:length(pcl_df$V1)] = 184:length(pcl_df$V1)
rownames(pcl_df) <- full_names
pcl_df <- pcl_df %>%
  select(-V1) %>%
  as.matrix() %>% t() %>%
  as_tibble()
pcl_df$participant_id <- pcl_df$study_id %>% str_sub(1,6)
# pcl_df %>% select(study_id, participant_id) %>% View()
# View(pcl_df[,c(1, 1859)])
# View(pcl_df[,1:100])

# get top 40 pathways
pathways <- grep(": ", names(full_names), perl = T)
pcl_pathways_df <- pcl_df %>% filter(omni_comparison == "comparison") %>%
  select(participant_id, collection_type, pathways)
pcl_pathways_df$participant_id[16] <- paste(pcl_pathways_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
pcl_pathways_df$participant_id[18] <- paste(pcl_pathways_df$participant_id[18], "b", sep = "_")
pcl_pathways_df <- pcl_pathways_df %>%
  gather("pathway", "val", -c(participant_id, collection_type))
pcl_pathways_df$val <- as.numeric(pcl_pathways_df$val)
pcl_pathways_df$val <- log2(pcl_pathways_df$val + 1) # allows for 0's
path <- pcl_pathways_df %>%
  mutate(pathway = full_names[pathway]) %>%
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
pdf(file.path(graph_dir, paste(desc, "Pathways_By_Collection_Diff.pdf", sep = "_")), width = 18, height = 12)
comp_pathways %>% ggplot(aes(x = reorder(pathway, -collection_diff), y = collection_diff)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance Diff (omni - reg)") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

# plot top pathways
pcl_pathways_df2 <- pcl_pathways_df %>%
  filter(pathway %in% top_12_path) %>%
  group_by(pathway, collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
# View(pcl_pathways_df)
pdf(file.path(graph_dir, paste(desc, "Top_40_Pathways_Side.pdf", sep = "_")), width = 18, height = 12)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(pathway, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pdf(file.path(graph_dir, paste(desc, "Top_40_Pathways_Stack.pdf", sep = "_")), width = 18, height = 12)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(pathway, -mean_val), y = mean_val, fill = collection_type)) +
  geom_bar(position = "stack", stat = "identity") +
  xlab("Pathway") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
pcl_pathways_df2 <- pcl_pathways_df %>%
  filter(pathway %in% top_12_path) %>%
  group_by(collection_type) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(mean_val)
pdf(file.path(graph_dir, paste(desc, "Pathways_Aggregated.pdf", sep = "_")), width = 6, height = 4)
pcl_pathways_df2 %>% ggplot(aes(x = reorder(collection_type, -mean_val), y = mean_val)) +
  geom_bar(position = "dodge", stat = "identity") +
  xlab("Collection Type") +
  ylab("Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()