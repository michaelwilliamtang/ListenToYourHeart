# 6/10/20
# Compares omni vs reg for batch metaphlan

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/Unbatched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
metaphlan_df <- met_taxa_df

omni_reg_comp <- function(met_lev, met_reg, met_lab, sing_lab, N, comp_N, is_batch) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  desc1 <- paste("Omnigene vs Regular, Top", comp_N, met_lab, "By Collection Difference", sep = " ")
  desc2 <- paste("Omnigene vs Regular, Top", N, met_lab, "By Relative Abundance", sep = " ")
  
  met_phyla_df <- metaphlan_df %>% 
    filter(omni_comparison == "comparison" &
           level == met_lev &
           batch == is_batch) %>%
    mutate(taxa = str_replace(taxa, met_reg, ""))
  phy <- met_phyla_df %>%
    group_by(taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$taxa[1:N]
  
  # plot most diff phyla, omni vs reg
  comp_phyla <- met_phyla_df %>%
    group_by(taxa, collection_type) %>%
    summarize(mean_val = mean(val)) %>%
    spread(collection_type, mean_val) %>%
    mutate(collection_diff = omni - regular,
           collection_diff_abs = abs(collection_diff)) %>%
    arrange(desc(collection_diff_abs))
  diff_N_phy <- comp_phyla$taxa[1:comp_N]
  comp_phyla <- comp_phyla %>% filter(taxa %in% diff_N_phy)
  fp <- file.path(graph_dir2, paste(met_lab, "By_Collection_Diff.pdf", sep = "_"))
  pdf(fp,  width = 18, height = 12)
  gg <- comp_phyla %>% ggplot(aes(x = reorder(taxa, -collection_diff), y = collection_diff)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab(sing_lab) +
    ylab("Relative Abundance Diff (omni - reg)") +
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
    labs(title = desc1)
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
  
  # plot top phyla
  met_phyla_df2 <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_")) %>%
    group_by(taxa, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(met_phyla_df)
  fp <- file.path(graph_dir2, paste("Top", N, met_lab, "Stack.pdf", sep = "_"))
  pdf(fp,  width = 18, height = 12)
  gg <- met_phyla_df2 %>% ggplot(aes(x = study_id, y = mean_val, fill = taxa)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Sample") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
    guides(fill = guide_legend(title = sing_lab)) +
    labs(title = desc2)
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
  met_phyla_df2 <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    group_by(collection_type, taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  fp <- file.path(graph_dir2, paste(met_lab, "Aggregated.pdf", sep = "_"))
  pdf(fp, width = 18, height = 12)
  gg <- met_phyla_df2 %>% ggplot(aes(x = collection_type, y = mean_val, fill = taxa)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Collection Type") +
    ylab("Relative Abundance") +
    theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
    guides(fill = guide_legend(title = sing_lab)) +
    labs(title = paste(desc2, ", Aggregated", sep = ""))
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot())
}

### comparing for metaphlan, phyla
omni_reg_comp(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, comp_N = 20, is_batch = F)

### comparing for metaphlan, genera
omni_reg_comp(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, comp_N = 20, is_batch = F)

### same, but batched
graph_dir <- "Graphs/Metaphlan/Batched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for metaphlan, phyla
omni_reg_comp(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, comp_N = 20, is_batch = T)

### comparing for metaphlan, genera
omni_reg_comp(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, comp_N = 20, is_batch = T)