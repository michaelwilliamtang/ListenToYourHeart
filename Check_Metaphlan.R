# 6/11/20
# Checks summation of unlogged metaphlan values and graphs distributions, per level

library(tidyverse)
data_dir <- "Data/Tidy"
summarize <- dplyr::summarize
ds1 <- "metaphlan"

check_distribution <- function(logged = F, met_lev, met_reg, met_lab, is_batch = F, coll_type = "omni") {
  print(met_lev)
  
  graph_dir2 <- "Graphs/Metaphlan/Check_Distribution"
  if (logged) graph_dir2 <- paste(graph_dir2, "Log", sep = "_")
  if (!dir.exists(graph_dir2)) dir.create(graph_dir2)
  load(file.path(data_dir, paste("Tidy_", ds1, ".RData", sep = "")))
  metaphlan_df <- met_taxa_df
  
  if (!is.na(is_batch)) metaphlan_df <- metaphlan_df %>%
    filter(batch == is_batch)
  if (!is.na(coll_type)) metaphlan_df <- metaphlan_df %>%
    filter(collection_type == coll_type)
  
  met_phyla_df <- metaphlan_df %>%
    filter(level == met_lev &
           omni_comparison == "comparison") %>%
    mutate(taxa = str_replace(taxa, met_reg, ""))
  met_phyla_df$batched <- "unbatched"
  met_phyla_df$batched[which(met_phyla_df$batch)] <- "batched"
  
  # unlog
  # met_phyla_df$val <- 2 ^ met_phyla_df$val - 1
  
  check_df <- met_phyla_df %>% 
    group_by(batched, collection_type, study_id) %>%
    summarize(total = sum(val)) %>%
    arrange(desc(total))
  write_tsv(check_df, file.path(graph_dir2, paste(met_lab, "Totals.txt", sep = "_")))
  
  check_df <- met_phyla_df %>% 
    group_by(batched, collection_type, study_id, taxa) %>%
    summarize(total = sum(val)) %>%
    arrange(desc(total))
  write_tsv(check_df, file.path(graph_dir2, paste(met_lab, "Taxa_Totals.txt", sep = "_")))
  
  pdf(file.path(graph_dir2, paste(met_lab, "Check_Distribution.pdf", sep = "_")), width = 18, height = 12)
  print(met_phyla_df %>% ggplot(aes(x = val, fill = taxa)) +
    geom_histogram(binwidth = 0.1, alpha = 0.5, position = "identity"))
  dev.off()
}


met_levs <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
met_regs <- c(".*k__", ".*p__", ".*c__", ".*o__", ".*f__", ".*g__", ".*s__")
met_labs <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
names(met_regs) <- met_levs
### comparing for metaphlan, phyla

for (i in 1:length(met_levs)) {
  check_distribution(logged = F, met_lev = met_levs[i], met_reg = met_regs[i], met_lab = met_labs[i],
                     is_batch = F, coll_type = "omni")
}