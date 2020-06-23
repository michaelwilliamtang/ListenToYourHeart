# 6/22/20
# Graph participant correlations from both iPOP1-2 and iPOP3 data, correlated against each other

library(tidyverse)

ncmr_cross_investigate <- function(source_dir) {
  print(source_dir)
  
  # get locs
  data_dir <- file.path("Data", source_dir)
  graph_dir <- file.path("Graphs", "RNA_Seq", paste(source_dir, "Comparison", sep = "_"))
  if (!dir.exists(graph_dir)) dir.create(graph_dir)
  
  # read data
  data <- read.table(file.path(data_dir, "ipop_all.txt"), stringsAsFactors = F)
  data <- data %>% select(-c(V2, V5))
  data$V1 <- gsub("_star_genome.vcf", "", data$V1)
  data$V3 <- gsub("_star_genome.vcf", "", data$V3)
  # filter specific
  unmatched <- c("SCGPM_Fiber-iPOP-1_C927Y_L1_unmatched_R1_SCGPM_Fiber-iPOP-1_C927Y_L1_unmatched_R2",
                 "SCGPM_Fiber-iPOP-2_C926V_L8_unmatched_R1_SCGPM_Fiber-iPOP-2_C926V_L8_unmatched_R2")
  data <- data %>% filter(!(V1 %in% unmatched) & !(V3 %in% unmatched))
  
  # read guides
  guide <- read.table(file.path(data_dir, "new_name_mapping_1-2.txt"), sep= "\t", header = 1, stringsAsFactors = F)
  guide2 <- read.table(file.path(data_dir, "new_name_mapping_3.txt"), sep= "\t", header = 1, stringsAsFactors = F)
  
  # add identifiers, merge
  guide$iPOP <- "1-2"
  guide2$iPOP <- "3"
  guide <- rbind(guide, guide2)
  
  id_index <- guide$sample.name %>% str_locate("_") - 1 # dynamically finds any length id
  guide$id <- guide$sample.name %>% substr(0, id_index)
  guide$id <- gsub("^69", "069", guide$id, perl = T) # regex
  
  # get missing
  missing_rows <- which(!(data$V1 %in% guide$new_full_name))
  write.csv(data$V[missing_rows], file.path(graph_dir, paste("Missing.csv", sep = "_")))
  
  # match data with guide
  samples <- guide$sample.name
  names(samples) <- guide$new_full_name
  ids <- guide$id
  names(ids) <- guide$new_full_name
  iPOP <- guide$iPOP
  names(iPOP) <- guide$new_full_name
  data <- data %>% filter(V1 %in% guide$new_full_name &
                            V3 %in% guide$new_full_name) %>%
    mutate(V2 = V3,
           corr = V4,
           id = ids[V1],
           sample = paste(iPOP[V1], samples[V1], sep = "_"), # also add identifers to samp names
           id2 = ids[V2],
           sample2 = paste(iPOP[V2], samples[V2], sep = "_")) %>%
    select(corr, id, sample, id2, sample2) # get rid of old cols
  
  # append reversed duplicate of the data
  data2 <- data[,]
  data2 <- data2 %>% mutate(tmp_id = id,
                            tmp_samp = sample,
                            id = id2,
                            sample = sample2,
                            id2 = tmp_id,
                            sample2 = tmp_samp) %>%
    select(-tmp_samp, -tmp_id)
  data <- rbind(data, data2)
  
  # write cleaned data
  write.csv(data, file.path(graph_dir, "Tidy_iPOP.csv"));
  
  # init cumu df
  mismatches_all_samps <- data.frame(samp = character(),
                                     total_corr = integer(),
                                     total_comp = integer(),
                                     mismatch = integer(),
                                     ratio = double(),
                                     low_corr = integer(),
                                     is_mismatch = logical())
  
  # plot per id
  all_ids <- data$id %>% unique()
  for (curr in all_ids) {
    print(curr)
    id_df <- data %>% filter(id == curr)
    gg <- id_df %>% ggplot(aes(x = sample, y = corr)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes(color = id2)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("correlation")
    ggsave(file.path(graph_dir, paste0("Correlations_", curr, ".pdf")), plot = last_plot())
    
    # init cumu df
    mismatches_by_samp <- data.frame(samp = character(),
                                     total_corr = integer(),
                                     total_comp = integer(),
                                     mismatch = integer(),
                                     ratio = double(),
                                     low_corr = integer(),
                                     is_mismatch = logical())
    
    # check per samp
    all_samples <- id_df$sample %>% unique()
    for (samp in all_samples) {
      # get ratio, categorized by corr strength and whether it's correlating an id with itself (non-mismatch)
      samp_df <- id_df %>% filter(sample == samp)
      high_corr <- samp_df %>% filter(corr > 0.65)
      low_corr <- samp_df %>% filter(corr <= 0.65)
      mismatch <- high_corr %>% filter(id != id2)
      ratio <- nrow(mismatch) / nrow(high_corr)
      
      if (ratio > 0.6 | is.na(ratio)) {
        is_mismatch <- T
      } else is_mismatch <- F
      
      # append to cumu per samp
      row <- data.frame(samp = samp,
                        total_corr = nrow(high_corr),
                        total_comp = nrow(samp_df),
                        mismatch = nrow(mismatch),
                        ratio = ratio,
                        low_corr = nrow(low_corr),
                        is_mismatch = is_mismatch)
      mismatches_by_samp <- mismatches_by_samp %>% rbind(row)
    }
    
    # write and append cumu per samp to cumu df
    write.csv(mismatches_by_samp, file.path(graph_dir, paste0("Mismatches_", curr, ".csv")))
    mismatches_all_samps <- mismatches_all_samps %>% rbind(mismatches_by_samp)
  }
  
  # write cumu df
  write.csv(mismatches_all_samps, file.path(graph_dir, "Mismatches_All.csv"))
  
  # build and write labeled guide
  guide2 <- guide
  guide2$is_mismatch <- mismatches_all_samps$is_mismatch[match(guide$sample.name, mismatches_all_samps$samp)]
  write.csv(guide2, file.path(graph_dir, "Guide_With_Mismatches.csv"))
  
  # plot by id
  gg <- data %>% ggplot(aes(x = id, y = corr)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = id2)) +
    ylab("correlation")
  ggsave(file.path(graph_dir, "Correlations_All.pdf"), plot = last_plot())
}

ncmr_cross_investigate("iPOP1-2-3")