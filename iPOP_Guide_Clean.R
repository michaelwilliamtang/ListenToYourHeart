# 6/22/20
# Cleans guide for iPOP3

data_dir <- file.path("Data", "iPOP3")

# read uncleaned guide
guide <- read.csv(file.path(data_dir, "iPOP3-guide.csv"), header = F)

# process id strings
IDs <- substring(as.character(guide[,1]), 2, 9)
IDs_last_8 <- substring(as.character(guide[,1]), 10)
# BiocManager::install("Biostrings")
library(Biostrings)
library(stringi)
dna = DNAStringSet(stri_reverse(IDs_last_8))
IDs_last_8_rc <- complement(dna)

# add/rename cols and write
file_names <- paste0(IDs, IDs_last_8_rc, "_R1_concatenated_", IDs, IDs_last_8_rc, "_R2_concatenated")
guide[,1] <- file_names
colnames(guide) <- c("new_full_name", "sample.name")
write.table(guide, file.path(data_dir, "new_name_mapping.txt"), sep = "\t", row.names = F)