# Reassign taxonomy for NWT seed microbiomes
# 1/31/2022 - Use latest SILVA and UNITE releases!
library(dada2)
library(dplyr)
library(tidyr)
seqs <- getSequences(system.file("extdata", "example_seqs.fa", package="dada2"))

#### ITS ####
seed_seqtab_ITS <- readRDS("Seeds/seqtab.nochimITS.rds")
ITS_seqs <- colnames(seed_seqtab_ITS)
ITS_seqs <- assignTaxonomy(ITS_seqs, 
                           "../../db_files/dada2/sh_general_release_dynamic_10.05.2021.fasta", 
                           tryRC = TRUE, 
                           multithread = TRUE)

# Flip table
seqtab.t <- as.data.frame(t(seed_seqtab_ITS))

# Pull out ASV repset
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
taxonomy <- as.data.frame(ITS_seqs)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
taxonomy_for_mctoolsr <- unite_(taxonomy, "taxonomy", 
                                c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"),
                                sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the taxonomy dataframe for the writeRepSetFasta function
taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file
writeRepSetFasta(taxonomy_for_fasta, "Seeds/ITSrepset.fasta")

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ASV <- NULL 
names(seqtab_wTax)[1] = "#ASV_ID"

# Set name of table in mctoolsr format and save
suppressWarnings(write.table(seqtab_wTax, "Seeds/seqtab_wTax_mctoolsr_ITS.txt", sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = "Seeds/seqtab_final_ITS.txt",
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(ITS_seqs, file = "Seeds/tax_final_ITS.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)



#### 16S ####
seed_seqtab_16S <- readRDS("Seeds/seqtab.nochim16S.rds")
seqs_16S <- colnames(seed_seqtab_16S)
seqs_16S <- assignTaxonomy(seqs_16S, 
                              "../../db_files/dada2/silva_nr99_v138.1_train_set.fa", 
                              tryRC = TRUE, 
                              multithread = TRUE)

# Flip table
seqtab.t <- as.data.frame(t(seed_seqtab_16S))

# Pull out ASV repset
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
taxonomy <- as.data.frame(seqs_16S)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
taxonomy_for_mctoolsr <- unite_(taxonomy, "taxonomy", 
                                c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"),
                                sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the taxonomy dataframe for the writeRepSetFasta function
taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ASV_ID"), 
        sep = ";", remove = FALSE) %>%
  unite("name", c("ASV_ID", "TaxString"), 
        sep = " ", remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file
writeRepSetFasta(taxonomy_for_fasta, "Seeds/16Srepset.fasta")

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ASV <- NULL 
names(seqtab_wTax)[1] = "#ASV_ID"

# Set name of table in mctoolsr format and save
suppressWarnings(write.table(seqtab_wTax, "Seeds/seqtab_wTax_mctoolsr_16S.txt", sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = "Seeds/seqtab_final_16S.txt",
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(seqs_16S, file = "Seeds/tax_final_16S.txt", 
            sep = "\t", row.names = TRUE, col.names = NA)
