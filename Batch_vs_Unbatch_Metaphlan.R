# 6/10/20
# Compares batch vs unbatch for omni or reg metaphlan

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/Omni_Batched_vs_Unbatched"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"

load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
metaphlan_df <- met_taxa_df

batch_unbatch_comp <- function(met_lev, met_reg, met_lab, sing_lab, N, comp_N, coll_type) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  met_phyla_df <- metaphlan_df %>% 
    filter(omni_comparison == "comparison" &
             level == met_lev &
             collection_type == coll_type) %>%
    mutate(taxa = str_replace(taxa, met_reg, ""))
  met_phyla_df$batched <- "unbatched"
  met_phyla_df$batched[which(met_phyla_df$batch)] <- "batched"
  phy <- met_phyla_df %>%
    group_by(taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$taxa[1:N]
  
  # plot most diff phyla, batch vs unbatch
  comp_phyla <- met_phyla_df %>%
    group_by(taxa, batched) %>%
    summarize(mean_val = mean(val)) %>%
    spread(batched, mean_val) %>%
    mutate(batch_diff = batched - unbatched,
           batch_diff_abs = abs(batch_diff)) %>%
    arrange(desc(batch_diff_abs))
  diff_N_phy <- comp_phyla$taxa[1:comp_N]
  comp_phyla <- comp_phyla %>% filter(taxa %in% diff_N_phy)
  pdf(file.path(graph_dir2, paste(met_lab, "By_Batch_Diff.pdf", sep = "_")), width = 18, height = 12)
  print(comp_phyla %>% ggplot(aes(x = reorder(taxa, -batch_diff), y = batch_diff)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab(sing_lab) +
          ylab("Relative Abundance Diff (batched - unbatched)") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)))
  dev.off()
  
  # plot top phyla
  met_phyla_df2 <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, batched, sep = "_")) %>%
    group_by(taxa, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(met_phyla_df)
  pdf(file.path(graph_dir2, paste("Top", N, met_lab, "Stack.pdf", sep = "_")), width = 18, height = 12)
  print(met_phyla_df2 %>% ggplot(aes(x = study_id, y = mean_val, fill = taxa)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Sample") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)))
  dev.off()
  met_phyla_df2 <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    group_by(batched, taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  pdf(file.path(graph_dir2, paste(met_lab, "Aggregated.pdf", sep = "_")), width = 18, height = 12)
  print(met_phyla_df2 %>% ggplot(aes(x = batched, y = mean_val, fill = taxa)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Collection Type") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)))
  dev.off()
}

### comparing for metaphlan, phyla
batch_unbatch_comp(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, 
              comp_N = 20, coll_type = "omni")

### comparing for metaphlan, genera
batch_unbatch_comp(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, 
              comp_N = 20, coll_type = "omni")

### now, same but with regular
graph_dir <- "Graphs/Metaphlan/Regular_Batched_vs_Unbatched"
if (!dir.exists(graph_dir)) dir.create(graph_dir)

### comparing for metaphlan, phyla
batch_unbatch_comp(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, 
              comp_N = 20, coll_type = "regular")

### comparing for metaphlan, genera
batch_unbatch_comp(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, 
              comp_N = 20, coll_type = "regular")