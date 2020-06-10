# 6/9/20
# Create NMDS graph for metaphlan

library(tidyverse)
library(goeveg)
library(vegan)
data_dir <- "Data"
graph_dir <- "Graphs/Metaphlan/Batched_Omni_vs_Regular"
summarize <- dplyr::summarize
desc <- "Phyla"
set.seed(10)

# library(icesTAF)
# library(metaMA)
# library(lme4)
# library(gtools)

# read metaphlan
metaphlan_df <- read.table(file.path(data_dir, "hf-metaphlan-final-combat.pcl"), sep = "\t", stringsAsFactors = FALSE, row.names = 1)
metaphlan_df <- metaphlan_df %>%
  as.matrix() %>% t() %>%
  as_tibble()

# get top 5 phyla
phyla <- grep(".*p__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_phyla_df <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(study_id, participant_id, collection_type, phyla)
colnames(met_phyla_df) <- met_phyla_df %>% colnames() %>%
  str_replace(".*p__", "")
# met_phyla_df$participant_id[16] <- paste(met_phyla_df$participant_id[16], "b", sep = "_") # anomaly participant (dup)
# met_phyla_df$participant_id[18] <- paste(met_phyla_df$participant_id[18], "b", sep = "_")
met_phyla_df <- met_phyla_df %>%
  gather("taxa", "val", -c(study_id, participant_id, collection_type))
met_phyla_df$val <- as.numeric(met_phyla_df$val)
phy <- met_phyla_df %>%
  group_by(taxa) %>%
  summarize(mean_val = mean(val)) %>%
  arrange(desc(mean_val))
top_5_phy <- phy$taxa[1:5]

# subset data and metadata
taxa <- grep(".*__[A-Za-z]+$", colnames(metaphlan_df), perl = T)
met_data <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(taxa)
met_data <- met_data %>% 
  sapply(as.numeric) %>%
  as.matrix()
met_metadata <- metaphlan_df %>% filter(omni_comparison == "comparison") %>%
  select(-taxa)

# get nmds
pdf(file.path(graph_dir, paste(desc, "Elbow_Plot.pdf", sep = "_")), width = 18, height = 12)
dimcheckMDS(met_data, distance = "bray", k = 10, trymax = 20,
            autotransform = TRUE)
dev.off()
nmds <- metaMDS(met_data, distance = "bray", k = 2, trymax = 500, na.rm = T, autotransform = F)

# envfit scores, write dimensions and vector wts
nmds_scores <- as.data.frame(scores(nmds))  
nmds_scores <- cbind(nmds_scores, met_metadata)
write.table(nmds_scores, row.names = F, file = file.path(graph_dir, "dimensions.txt"), sep = "\t", quote = FALSE)
envf <- envfit(nmds, met_data, perm = 999, na.rm = T)
vector.weights <-  scores(envf, display = "vectors")
padjusted <- p.adjust(envf$vectors$pvals, method = "fdr")
vector.weights <- cbind(analyte = rownames(vector.weights), vector.weights, padjust = padjusted)
write.table(vector.weights, row.names = F, file = file.path(graph_dir, "vector_weights.txt"), sep = "\t", quote = FALSE)

gg <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, fill = participant_id, color = participant_id, shape = collection_type)) + 
  geom_point(size=3) +  # Set color to vary based on subject
  labs(title = desc)
plot(gg)
ggsave(file.path(graph_dir, paste(desc, "NMDS.pdf")), plot = last_plot())


# gg <- ggplot(data=nmds_scores, aes(x = NMDS1, y = NMDS3, fill = study_id, color = study_id, shape = collection_type)) + 
#   geom_point(size = 3) +  # Set color to vary based on subject
#   #labs(title="Metagenome_phylum", subtitle="All fiber types", y="NMDS2", x="NMDS1") +
#   theme()#legend.position="none"
# plot(gg)
# ggsave(file.path(graph_dir, paste(desc, "NMDS1-3.pdf")), plot = last_plot())
# 
# 
# gg <- ggplot(data=nmds_scores, aes(x = NMDS2,y = NMDS3, fill = study_id, color = study_id, shape = collection_type)) + 
#   geom_point(size = 3) +  # Set color to vary based on subject
#   #labs(title="Metagenome_phylum", subtitle="All fiber types", y="NMDS2", x="NMDS1") +
#   theme()#legend.position="none"
# plot(gg)
# ggsave(file.path(graph_dir, paste(desc, "NMDS2-3.pdf")), plot = last_plot())
