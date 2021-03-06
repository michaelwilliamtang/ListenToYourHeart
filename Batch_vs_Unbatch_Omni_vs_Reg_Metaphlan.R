# 6/11/20
# Compares batch omni vs unbatch omni vs batch reg vs omni reg metaphlan

library(tidyverse)
data_dir <- "Data/Tidy"
graph_dir <- "Graphs/Metaphlan/Batched_vs_Unbatched_Omni_vs_Regular_Scaled"
if (!dir.exists(graph_dir)) dir.create(graph_dir)
summarize <- dplyr::summarize
ds1 <- "metaphlan"

load(file.path(data_dir, paste("Tidy_Scaled_", ds1, ".RData", sep = "")))
metaphlan_df <- met_taxa_df

batch_unbatch_omni_reg_comp <- function(met_lev, met_reg, met_lab, sing_lab, N, comp_N) {
  graph_dir2 <- file.path(graph_dir, paste("Top", N, met_lab, sep = "_"))
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  
  # desc1 <- paste("Omnigene vs Regular and Batch Corrected vs Uncorrected\nTop", 
  #   comp_N, met_lab, "By Collection Difference", sep = " ")
  desc2 <- paste("Omnigene vs Regular and Batch Corrected vs Uncorrected\nTop", 
                 N, met_lab, "By Relative Abundance", sep = " ")
  
  met_phyla_df <- metaphlan_df %>% 
    filter(omni_comparison == "comparison" &
             level == met_lev) %>%
    mutate(taxa = str_replace(taxa, met_reg, ""))
  met_phyla_df$batched <- "unbatched"
  met_phyla_df$batched[which(met_phyla_df$batch)] <- "batched"
  phy <- met_phyla_df %>%
    group_by(taxa) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(desc(mean_val))
  top_N_phy <- phy$taxa[1:N]
  
  # plot top phyla
  met_phyla_df2 <- met_phyla_df %>%
    filter(taxa %in% top_N_phy) %>%
    mutate(study_id = paste(study_id, collection_type, sep = "_")) %>%
    mutate(study_id = paste(study_id, batched, sep = "_")) %>%
    group_by(taxa, study_id) %>%
    summarize(mean_val = mean(val)) %>%
    arrange(mean_val)
  # View(met_phyla_df)
  fp <- file.path(graph_dir2, paste("Top", N, met_lab, "Stack.pdf", sep = "_"))
  pdf(fp, width = 18, height = 12)
  gg <- met_phyla_df2 %>% ggplot(aes(x = study_id, y = mean_val, fill = taxa)) +
          geom_bar(position = "stack", stat = "identity") +
          xlab("Sample") +
          ylab("Relative Abundance") +
          theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
          guides(fill = guide_legend(title = sing_lab)) +
          labs(title = desc2)
  plot(gg)
  dev.off()
  # ggsave(fp, plot = last_plot(), width = 18, height = 12)
  # met_phyla_df2 <- met_phyla_df %>%
  #   filter(taxa %in% top_N_phy) %>%
  #   group_by(batched, taxa) %>%
  #   summarize(mean_val = mean(val)) %>%
  #   arrange(mean_val)
  # pdf(file.path(graph_dir2, paste(met_lab, "Aggregated.pdf", sep = "_")), width = 18, height = 12)
  # print(met_phyla_df2 %>% ggplot(aes(x = batched, y = mean_val, fill = taxa)) +
  #         geom_bar(position = "stack", stat = "identity") +
  #         xlab("Collection Type") +
  #         ylab("Relative Abundance") +
  #         theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
  #         guides(fill = guide_legend(title = sing_lab)))
  # dev.off()
}

### comparing for metaphlan, phyla
batch_unbatch_omni_reg_comp(met_lev = "phylum", met_reg = ".*p__", met_lab = "Phyla", sing_lab = "Phylum", N = 5, 
                   comp_N = 20)

### comparing for metaphlan, genera
batch_unbatch_omni_reg_comp(met_lev = "genus", met_reg = ".*g__", met_lab = "Genera", sing_lab = "Genus", N = 20, 
                   comp_N = 20)