# Methylotrophic, hydrogenotrophhic, acetotrophic methanogenesis genes comparative analysis
# Compare the methanogenesis habitats from Lyu 2018
# For review paper
# by Cliff Bueno de Mesquita, JGI, Spring 2022
# Revised October 2022 for revisions for MMBR

#### Retrieving Data ####
# To retrieve data from IMG/M, search all fields for the environment names
# Filter out SPAdes duplicates, metatranscriptomes, Bacteria, Archaea
# To retrieve functional data (KO), use the statistical analysis tools.
# ^ this requires you to split up the genome set into at least 2 batches, but you can still get the raw data in the output!

#### Setup ####
suppressWarnings(suppressMessages(library(readxl))) # For read_xlsx
suppressWarnings(suppressMessages(library(janitor))) # For cleaning
suppressWarnings(suppressMessages(library(cowplot))) # For multipanel
suppressWarnings(suppressMessages(library(plyr))) # For data manipulation
suppressWarnings(suppressMessages(library(tidyverse))) # For data manipulation
suppressWarnings(suppressMessages(library(reshape2))) # For melting
suppressWarnings(suppressMessages(library(vegan))) # For analysis
suppressWarnings(suppressMessages(library(car))) # For leveneTest
suppressWarnings(suppressMessages(library(PMCMRplus))) # For Nemenyi posthoc test
suppressWarnings(suppressMessages(library(indicspecies))) # For multipatt
suppressWarnings(suppressMessages(library(scales))) # For muted
suppressWarnings(suppressMessages(library(DESeq2))) # For normalization
suppressWarnings(suppressMessages(library(FSA))) # For se
suppressWarnings(suppressMessages(library(mctoolsr))) # For taxonomic analysis
suppressWarnings(suppressMessages(library(cowplot))) # For multipanel
suppressWarnings(suppressMessages(library(plotly))) # For interactive graphs
suppressWarnings(suppressMessages(library(RColorBrewer))) # For color palettes
suppressWarnings(suppressMessages(library(dendextend))) # For dendrogram plots
suppressWarnings(suppressMessages(library(viridis))) # For viridis palette
suppressWarnings(suppressMessages(library(gplots))) # For heatmaps
suppressWarnings(suppressMessages(library(maps))) # For geographic maps
suppressWarnings(suppressMessages(library(mapproj))) # For geographic maps
suppressWarnings(suppressMessages(library(magrittr))) # For setting column names
suppressWarnings(suppressMessages(library(writexl))) # For writing Excel file
suppressWarnings(suppressMessages(library(plotrix))) # For standard error
suppressWarnings(suppressMessages(library(ggh4x))) # For nested facet wrapping

setwd("~/Desktop/Review")
options(max.print = 20000000)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullj <- function(df) df[chull(df$Axis01j, df$Axis02j),]
`%notin%` <- Negate(`%in%`)



#### Initial Setup ####
# Original data input and processing. Don't need to repeat. Can start at "START HERE"

# Functional table (KO)
f <- read.delim("lyu_ko_statanalysis_09-feb-2022/stats_input.txt") %>%
  filter(UniqueID != "GroupID") %>%
  separate(UniqueID, into = c("Text", "KO"), sep = ":") %>%
  select(-Text) %>%
  column_to_rownames(var = "KO") %>%
  replace(is.na(.), 0) %>%
  mutate_if(is.character, as.numeric)

f_comm <- as.data.frame(t(f)) %>%
  arrange(row.names(.))

# Metadata
f_meta <- read.delim("lyu_habitats_metadata.txt") %>%
  arrange(taxon_oid) %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>%
  mutate_if(is.character, as.factor)
sum(f_meta$sampleID != rownames(f_comm))
f_meta$richness_KO = specnumber(f_comm)
sum(f_meta$richness_KO == 0) # 19 zeroes, remove.

f_comm <- f_comm %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(richness_KO = f_meta$richness_KO) %>%
  filter(richness_KO > 0) %>%
  select(-richness_KO)
f_meta <- subset(f_meta, richness_KO > 0)
sum(f_meta$sampleID != f_comm$sampleID)
rownames(f_comm) <- f_comm$sampleID
f_comm$sampleID <- NULL

# Make ecosystem variable
f_meta$Ecosystem <- NA
for (i in 1:nrow(f_meta)) {
  if(length(grep("landfill", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("landfill", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Landfill"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("sewage", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("sewage", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Sewage treatment"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("wetland", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("wetland", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Wetland (freshwater)"
  }
}

for (i in 1:nrow(f_meta)) {
  if(f_meta$Ecosystem[i] == "Wetland (freshwater)" &
     length(grep("coastal", f_meta$Study.Name[i]) > 0) |
     length(grep("Coastal", f_meta$Study.Name[i]) > 0) |
     length(grep("coastal", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0) |
     length(grep("R2A", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0) |
     length(grep("SF2", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0) |
     length(grep("Natural and restored wetland microbial communities from the San Francisco Bay, California, USA, that impact long-term carbon sequestration", f_meta$Study.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Wetland (coastal)"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("rice", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("rice", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Rice field"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("ocean", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("ocean", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Ocean"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("human gut", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("human gut", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Human gut"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("ruminant", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("ruminant", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Sheep gut"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("cow", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("cow", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Cow gut"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("termite", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("termite", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Termite gut"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("hydrothermal", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("hydrothermal", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Hydrothermal vent"
  }
}

for (i in 1:nrow(f_meta)) {
  if(length(grep("hypersaline", f_meta$Study.Name[i], ignore.case = T) > 0) |
     length(grep("hypersaline", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0) |
     length(grep("R1_", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0) |
     length(grep("R2_", f_meta$Genome.Name...Sample.Name[i], ignore.case = T) > 0)) {
    f_meta$Ecosystem[i] <- "Hypersaline"
  }
}

f_meta <- f_meta %>%
  mutate(Ecosystem = as.factor(Ecosystem),
         Hypothesis = recode_factor(Ecosystem,
                                    "Cow gut" = "H2",
                                    "Human gut" = "H2",
                                    "Hydrothermal vent" = "H2",
                                    "Hypersaline" = "Me",
                                    "Landfill" = "H2, Ac",
                                    "Ocean" = "Me",
                                    "Rice field" = "H2, Ac",
                                    "Sewage treatment" = "H2, Ac",
                                    "Sheep gut" = "H2",
                                    "Termite gut" = "H2",
                                    "Wetland (coastal)" = "Me",
                                    "Wetland (freshwater)" = "H2, Ac"))

table(f_meta$Ecosystem)
table(f_meta$Hypothesis)

# Sample size counts and plot
system <- as.data.frame(table(f_meta$Ecosystem))
ggplot(system, aes(reorder(Var1, Freq, mean), Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Ecosystem", 
       y = "Count") +
  coord_flip() +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))

# Sample map
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", color = "white", size = 0.5) +
  coord_map(projection = "mollweide") +
  labs(x = "Longitude",
       y = "Latitude") + 
  theme_map() +
  geom_point(data = f_meta, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5)


# KO richness by sequencing depth
ggplot(f_meta, aes(Genome.Size.....assembled, richness_KO)) +
  geom_point(size = 1.5, alpha = 0.25) +
  geom_smooth() +
  labs(x = "Assembled genome size", 
       y = "# KOs") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# DESeq Normalization
dds <- DESeqDataSetFromMatrix(countData = t(f_comm) + 1,
                              colData = f_meta,
                              design = ~ 1)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
f_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))

# Save files (then you don't have to repeat this part)
saveRDS(f_comm_DESeq, "f_comm_DESeq.rds")
saveRDS(f_comm, "f_comm.rds")
saveRDS(f_meta, "f_meta.rds")


#### START HERE ####
#### Functional  ####
# Load tables
f_comm <- readRDS("f_comm.rds")
f_comm_DESeq <- readRDS("f_comm_DESeq.rds")
f_meta <- readRDS("f_meta.rds") %>%
  mutate(sampleID = paste("X", IMG.Genome.ID, sep = ""))
write_xlsx(f_meta, 
          "~/Desktop/Review/habitats_permissions.xlsx",
           format_headers = F)

f_meta_lyu <- f_meta %>%
  filter(Ecosystem != "Hypersaline") %>%
  filter(Ecosystem != "Wetland (coastal)")
# write_xlsx(f_meta_lyu, 
#           "~/Desktop/Review/lyu_habitats_permissions.xlsx",
#           format_headers = F)
  

# Check IDs
sum(f_meta$sampleID != rownames(f_comm_DESeq))

#### _Differential Abundance ####
# Differential abundance analaysis on all genes, with DESeq2
dds_input <- DESeqDataSetFromMatrix(countData = t(f_comm) + 1,
                                    colData = f_meta,
                                    design = ~ Ecosystem)
dds <- DESeq(dds_input)
res <- results(dds)
res
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Ecosystem_Hydrothermal.vent_vs_Cow.gut", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

# p = 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2))

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-11,12))
title("317 Genes up- or down-regulated in metatranscriptomes")
dev.off()
plotCounts(dds, gene=which.min(res$padj), intgroup="Type")
res_df <- data.frame("KO" = colnames(f_comm),
                     "log2FoldChange" = res$log2FoldChange,
                     "se" = res$lfcSE,
                     "Padj" = res$padj) %>%
  filter(Padj < 0.05)

# Get list of KOs for MG and MT
MT_genes <- filter(res_df, log2FoldChange > 0)
MG_genes <- filter(res_df, log2FoldChange < 0)
# write.csv(MT_genes, "MT_Upreg_Genes.csv")



#### _Methylotrophic Genes ####

# Plot genes by ecosystem
# mtgB not present
gene_plot <- data.frame("Ecosystem" = f_meta$Ecosystem,
                        "Hypothesis" = f_meta$Hypothesis,
                        "log2DESeqmttB" = log2(f_comm_DESeq$K14083),
                        "log2DESeqmttC" = log2(f_comm_DESeq$K14084),
                        "log2DESeqmtbB" = log2(f_comm_DESeq$K16178),
                        "log2DESeqmtbC" = log2(f_comm_DESeq$K16179),
                        "log2DESeqmtmB" = log2(f_comm_DESeq$K16176),
                        "log2DESeqmtmC" = log2(f_comm_DESeq$K16177),
                        "log2DESeqmtbA" = log2(f_comm_DESeq$K14082),
                        "log2DESeqmtaB" = log2(f_comm_DESeq$K04480),
                        "log2DESeqmtaC" = log2(f_comm_DESeq$K14081),
                        "log2DESeqmtaA" = log2(f_comm_DESeq$K14080),
                        "log2DESeqmtsA" = log2(f_comm_DESeq$K16954),
                        "Database" = "IMG")

leveneTest(log2DESeqmttB ~ Ecosystem, data = gene_plot)
m <- aov(log2DESeqmttB ~ Ecosystem, data = gene_plot)
summary(m)
shapiro.test(m$residuals)
TukeyHSD(m)
kruskal.test(log2DESeqmttB ~ Ecosystem, data = gene_plot)
kwAllPairsNemenyiTest(log2DESeqmttB ~ Ecosystem, data = gene_plot)

ggplot(gene_plot, aes(reorder(Hypothesis, log2DESeqmttB, mean), log2DESeqmttB)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Substrate", 
       y = "mttB abundance (log2 DESeq count)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))

ggplot(gene_plot, aes(reorder(Hypothesis, log2DESeqmttC, mean), log2DESeqmttC)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Substrate", 
       y = "mttC abundance (log2 DESeq count)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))



#### _Multipanel Figure ####
gene_plot_long <- gene_plot %>%
  pivot_longer(c("log2DESeqmttB", "log2DESeqmttC", "log2DESeqmtbB", "log2DESeqmtbC",
                 "log2DESeqmtmB", "log2DESeqmtmC", "log2DESeqmtbA", "log2DESeqmtaB",
                 "log2DESeqmtaC", "log2DESeqmtaA", "log2DESeqmtsA"), 
               names_to = "Gene", values_to = "Abundance")

pdf("~/Desktop/Review/Gene_Ecosystem_Counts.pdf", width = 10, height = 8)
ggplot(gene_plot_long, aes(reorder(Ecosystem, Abundance, mean), Abundance)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "Log2 DESeq Abundance") +
  facet_wrap(~ Gene, ncol = 3) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))
dev.off()

# Me, H2, Ac, select one gene for each
# Also include mcrA, K00399
# Ac: cdhD, K00194
# H2: mtrA, K00577
# M: mtaA, K14080
# TMA: mttC, K14084
# DMA: mtbC, K16179
# MMA: mtmC, K16177
# DMS/MeSH/MMPA: mtsA, K16954
gene_plot2 <- data.frame("Ecosystem" = f_meta$Ecosystem,
                        "Hypothesis" = f_meta$Hypothesis,
                        "log2SESeqmcrA" = log2(f_comm_DESeq$K00399),
                        "log2DESeqcdhD" = log2(f_comm_DESeq$K00194),
                        "log2DESeqmtrA" = log2(f_comm_DESeq$K00577),
                        "log2DESeqmttC" = log2(f_comm_DESeq$K14084),
                        "log2DESeqmtbC" = log2(f_comm_DESeq$K16179),
                        "log2DESeqmtmC" = log2(f_comm_DESeq$K16177),
                        "log2DESeqmtaA" = log2(f_comm_DESeq$K14080),
                        "log2DESeqmtsA" = log2(f_comm_DESeq$K16954))
gene_plot2_long <- gene_plot2 %>%
  pivot_longer(c("log2SESeqmcrA", "log2DESeqcdhD", "log2DESeqmtrA", "log2DESeqmttC", 
                 "log2DESeqmtbC", "log2DESeqmtmC", "log2DESeqmtaA", "log2DESeqmtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "log2SESeqmcrA" = "mcrA",
                              "log2DESeqcdhD" = "cdhD",
                              "log2DESeqmtrA" = "mtrA",
                              "log2DESeqmttC" = "mttC",
                              "log2DESeqmtbC" = "mtbC", 
                              "log2DESeqmtmC" = "mtmC", 
                              "log2DESeqmtaA" = "mtaA",
                              "log2DESeqmtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "mtrA" = "H2/CO2 (mtrA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)"))

pdf("~/Desktop/Review/Gene_Ecosystem_Counts.pdf", width = 10, height = 8)
ggplot(gene_plot2_long, aes(Ecosystem, Abundance, colour = Pathway)) +
  geom_boxplot(outlier.size = 1) +
  labs(x = "Ecosystem", 
       y = "Log2 DESeq Abundance") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))
dev.off()

lyu_only <- gene_plot2_long %>%
  filter(Ecosystem != "Hypersaline") %>%
  filter(Ecosystem != "Wetland (coastal)") %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

ggplot(lyu_only, aes(Ecosystem, Abundance, colour = Pathway)) +
  geom_boxplot(outlier.size = 1) +
  labs(x = "Ecosystem", 
       y = "Log2 DESeq Abundance") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))

lyu_only_summary <- lyu_only %>%
  group_by(Hypothesis, Ecosystem, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

pdf("Figure3_DESeq.pdf", width = 8, height = 4)
ggplot(lyu_only_summary, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                    group = Pathway),
                position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = expression(bold('Log'[2]~'DESeq Abundance'))) +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8))
dev.off()

png("Figure3_DESeq.png", width = 8, height = 4, units = "in", res = 300)
ggplot(lyu_only_summary, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = expression(bold('Log'[2]~'DESeq Abundance'))) +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8))
dev.off()



#### _MuSICC #####
for_musicc <- f_comm %>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(var = "KO")
# write.csv(for_musicc, "ko_abundances.csv", row.names = FALSE)

# This giving NaN errors
# Filter out samples with 0 USCGs!
uscg <- read_excel("13059_2015_610_MOESM1_ESM.xlsx", sheet = 2)
for_musicc_uscg <- for_musicc %>%
  filter(KO %in% uscg$KO)
with_uscg <- data.frame("uscg_count" = colSums(for_musicc_uscg[,2:2812])) %>%
  filter(uscg_count > 0)
for_musicc_filt <- for_musicc %>%
  select(KO, rownames(with_uscg))
# write.csv(for_musicc_filt, "ko_abundances_filt.csv", row.names = FALSE)

# Input
f_comm_MUSiCC <- read.csv("MUSiCC_normalized.csv") %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>%
  na.omit()

f_meta_MUSiCC <- f_meta %>%
  filter(sampleID %in% rownames(f_comm_MUSiCC)) %>%
  filter(Ecosystem != "Hypersaline") %>%
  filter(Ecosystem != "Wetland (coastal)")

f_comm_MUSiCC <- f_comm_MUSiCC %>%
  filter(rownames(.) %in% f_meta_MUSiCC$sampleID)

# Check IDs
sum(f_meta_MUSiCC$sampleID != rownames(f_comm_MUSiCC))

# Data frame
gene_plot2_m_lyu <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                               "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                               "MUSiCCmcrA" = f_comm_MUSiCC$K00399,
                               "MUSiCCcdhD" = f_comm_MUSiCC$K00194,
                               "MUSiCCmtrA" = f_comm_MUSiCC$K00577,
                               "MUSiCCmttC" = f_comm_MUSiCC$K14084,
                               "MUSiCCmtbC" = f_comm_MUSiCC$K16179,
                               "MUSiCCmtmC" = f_comm_MUSiCC$K16177,
                               "MUSiCCmtaA" = f_comm_MUSiCC$K14080,
                               "MUSiCCmtsA" = f_comm_MUSiCC$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent"))

#### __ (a) Stats ####
# Run a loop 
kruskal_results_genes <- as.data.frame(matrix(data = NA, 10, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 3:10) {
  k <- kruskal.test(gene_plot2_m_lyu[[i]] ~ gene_plot2_m_lyu$Ecosystem)
  kruskal_results_genes$Gene[i] <- names(gene_plot2_m_lyu)[i]
  kruskal_results_genes$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes$P[i] <- k$p.value
}
kruskal_results_genes <- kruskal_results_genes %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))

# Run manually (not needed anymore)
leveneTest(MUSiCCmcrA ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmcrA ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmcrA ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCcdhD ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCcdhD ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCcdhD ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmtrA ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmtrA ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmtrA ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmtaA ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmtaA ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmtaA ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmttC ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmttC ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmttC ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmtbC ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmtbC ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmtbC ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmtmC ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmtmC ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmtmC ~ Ecosystem, data = gene_plot2_m_lyu)

leveneTest(MUSiCCmtsA ~ Ecosystem, data = gene_plot2_m_lyu)
m <- aov(MUSiCCmtsA ~ Ecosystem, data = gene_plot2_m_lyu)
summary(m)
shapiro.test(m$residuals)
kruskal.test(MUSiCCmtsA ~ Ecosystem, data = gene_plot2_m_lyu)



#### __ (b) Graph ####
gene_plot2_long_m <- gene_plot2_m %>%
  pivot_longer(c("log2SESeqmcrA", "MUSiCCcdhD", "MUSiCCmtrA", "MUSiCCmttC", 
                 "MUSiCCmtbC", "MUSiCCmtmC", "MUSiCCmtaA", "MUSiCCmtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "log2SESeqmcrA" = "mcrA",
                              "MUSiCCcdhD" = "cdhD",
                              "MUSiCCmtrA" = "mtrA",
                              "MUSiCCmttC" = "mttC",
                              "MUSiCCmtbC" = "mtbC", 
                              "MUSiCCmtmC" = "mtmC", 
                              "MUSiCCmtaA" = "mtaA",
                              "MUSiCCmtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "mtrA" = "H2/CO2 (mtrA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)"))

lyu_only_m <- gene_plot2_long_m %>%
  filter(Ecosystem != "Hypersaline") %>%
  filter(Ecosystem != "Wetland (coastal)") %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

lyu_only_summary_m <- lyu_only_m %>%
  group_by(Hypothesis, Ecosystem, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

pdf("Figure3_MUSiCC.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (MUSiCC normalized)") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()

png("Figure3_MUSiCC.png", width = 9, height = 4, unit = "in", res = 300)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (MUSiCC normalized)") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()



#### Taxonomic ####
# Metadata ("mapping file") downloaded from IMG

# Tax table for mctoolsr
t <- read.delim("lyu_family_24-feb-2022/stats_input") %>%
  filter(UniqueID != "GroupID") %>%
  dplyr::rename(taxonomy = UniqueID) %>%
  mutate(ASV_ID = taxonomy) %>%
  select(ASV_ID, 3:ncol(.)-1, taxonomy)
names(t) <- abbreviate(names(t), minlength = 11)
table.fp <- "~/Desktop/Review/"
out_fp <- paste0(table.fp, "/IMG_families_mctoolsr.txt")
names(t)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(t, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Import
tax_table_fp <- file.path("~/Desktop/Review/IMG_families_mctoolsr.txt")
map_fp <- file.path("~/Desktop/Review/lyu_habitats_metadata.txt")
input = load_taxa_table(tax_table_fp, map_fp)
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""))

# Match the functional analysis (removed 19 zeroes)
input = filter_data(input,
                    filter_cat = "sampleID",
                    keep_vals = f_meta_lyu$sampleID)
input$map_loaded <- input$map_loaded %>%
  left_join(., f_meta_lyu[,25:28], by = "sampleID") %>%
  column_to_rownames(var = "sampleID") %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""))

# Check sequencing depth 
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded))
se(colSums(input$data_loaded))
input$map_loaded$count <- colSums(input$data_loaded)
ggplot(input$map_loaded, aes(reorder(Ecosystem, count, mean), count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "# Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

ggplot(input$map_loaded, aes(`Genome Size   * assembled`, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  labs(x = "Assembled genome size", 
       y = "Assigned family reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# For all data get rid of samples with less than 1000 family reads (n = 73)
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` < 1000)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))
mean(colSums(input$data_loaded))
se(colSums(input$data_loaded))
min(colSums(input$data_loaded))
max(colSums(input$data_loaded))

# Rename size column
input$map_loaded$GenomeSize = input$map_loaded$`Genome Size   * assembled`

# Phyla
tax_sum_phyla <- summarize_taxonomy(input, level = 2, report_higher_tax = T, relative = FALSE)
# Note there are no Bathyarchaeota or Verstraetearchaeota

# Summarize by family, extract methanogens, calculate CPM
tax_sum_family_wTax <- summarize_taxonomy(input, level = 5, report_higher_tax = T, relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
methano_wTax <- methano_wTax[!grepl("Plasmid", rownames(methano_wTax)),]
# Note this also contains a methanotroph Methermicoccaceae and a plasmid Methanobacteriales

# Without higher tax and no plasmid
tax_sum_family <- summarize_taxonomy(input, level = 5, report_higher_tax = FALSE, relative = FALSE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
methano <- subset(methano, rownames(methano) != "Methanobacteriales")
tax_sum_family_df <- t(methano) %>%
  as.data.frame() %>%
  mutate(`Candidatus Methanoperedenaceae` = (`Candidatus Methanoperedenaceae`*1000000)/input$map_loaded$GenomeSize,
         Methanobacteriaceae = (Methanobacteriaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocaldococcaceae = (Methanocaldococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocellaceae = (Methanocellaceae*1000000)/input$map_loaded$GenomeSize,
         Methanococcaceae = (Methanococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocorpusculaceae = (Methanocorpusculaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomassiliicoccaceae = (Methanomassiliicoccaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomicrobiaceae = (Methanomicrobiaceae*1000000)/input$map_loaded$GenomeSize,
         Methanonatronarchaeaceae = (Methanonatronarchaeaceae*1000000)/input$map_loaded$GenomeSize,
         Methanopyraceae = (Methanopyraceae*1000000)/input$map_loaded$GenomeSize,
         Methanoregulaceae = (Methanoregulaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosaetaceae = (Methanosaetaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosarcinaceae = (Methanosarcinaceae*1000000)/input$map_loaded$GenomeSize,
         Methanospirillaceae = (Methanospirillaceae*1000000)/input$map_loaded$GenomeSize,
         Methanothermaceae = (Methanothermaceae*1000000)/input$map_loaded$GenomeSize,
         Methanotrichaceae = (Methanotrichaceae*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanomicrobiales` = (`unclassified Methanomicrobiales`*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanosarcinales` = (`unclassified Methanosarcinales`*1000000)/input$map_loaded$GenomeSize) 

input$map_loaded <- input$map_loaded %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))
sum(input$map$sampleID != rownames(tax_sum_family_df))
input$map_loaded <- cbind(input$map_loaded, tax_sum_family_df)



#### (a) Stats ####
# Could do manually or put df back to mctoolsr format for built in function
# mctoolsr
tax_sum_family_methano_CPM <- t(tax_sum_family_df) %>%
  as.data.frame()

kw <- taxa_summary_by_sample_type(tax_sum_family_methano_CPM,
                                  metadata_map = input$map_loaded,
                                  type_header = "Ecosystem",
                                  test_type = "KW")

# This gives P but not X2, so run manually with for loop
kruskal_results_taxa <- as.data.frame(matrix(data = NA, 49, 3)) %>%
  set_names(c("Family", "X2", "P"))
for (i in 31:48) {
  k <- kruskal.test(input$map_loaded[[i]] ~ input$map_loaded$Ecosystem)
  kruskal_results_taxa$Family[i] <- names(input$map_loaded)[i]
  kruskal_results_taxa$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_taxa$P[i] <- k$p.value
}
kruskal_results_taxa <- kruskal_results_taxa %>%
  filter(is.na(Family) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))


#### (b) Graph ####
family_plot_long <- pivot_longer(input$map_loaded,
                                 cols = names(input$map_loaded)[30:48],
                                 names_to = "Family",
                                 values_to = "Abundance")
bars <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  filter(Family != "unclassified Methanomicrobiales",
         Family != "unclassified Methanosarcinales",
         Family != "Candidatus Methanoperedenaceae")

filter_bars <- bars %>%
  group_by(Family) %>%
  summarise(max_abund = max(mean_abund)) %>%
  filter(max_abund > 1)

bars_top <- bars %>%
  filter(Family %in% filter_bars$Family)

sort_bars <- bars_top %>%
  group_by(Family) %>%
  summarise(sum_abund = sum(mean_abund)) %>%
  arrange(sum_abund)

bars_top_sorted <- bars_top %>%
  mutate(Family = factor(Family,
                         levels = sort_bars$Family))

n_families <- bars %>%
  group_by(Hypothesis, Ecosystem) %>%
  summarise(n_fam = sum(mean_abund > 0),
            tot = sum(mean_abund)) %>%
  mutate(y = tot + 1)

nb.cols <- nrow(filter_bars)
mycolors <- colorRampPalette(brewer.pal(11, "Paired"))(nb.cols)

pdf("Figure4.pdf", width = 8.5, height = 4.5)
ggplot(bars_top_sorted, aes(Ecosystem, mean_abund, fill = Family)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.box.margin=margin(0,-5,0,15))
dev.off()

png("Figure4.png", width = 8.5, height = 4.5, unit = "in", res = 300)
ggplot(bars_top_sorted, aes(Ecosystem, mean_abund, fill = Family)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free") +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.box.margin=margin(0,-5,0,15))
dev.off()

ggplotly(g1)

#### (c) PCoA ####
# For just methanogen, so need to filter some more samples
input$data_loaded <- tax_sum_family_methano_CPM
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` == 0)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))
bc <- calc_dm(input$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input$map_loaded) - 1, eig = T)
eigenvals(pcoa)/sum(eigenvals(pcoa)) # 35.4, 14.5% variation explained
input$map_loaded$Axis01 <- scores(pcoa)[,1]
input$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input$map_loaded, c("Ecosystem"), find_hull)
ggplot(input$map_loaded, aes(Axis01, Axis02, colour = Ecosystem)) +
  geom_polygon(data = micro.hulls, aes(colour = Ecosystem, fill = Ecosystem),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "PC1: 35.4% Variation Explained", 
       y = "PC2: 14.5% Variation Explained",
       colour = "Habitat") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Very messy, probably not good as many samples but few taxa. Don't need to show in paper

# PERMANOVA
set.seed(1150)
input$map_loaded$Ecosystem <- as.factor(input$map_loaded$Ecosystem)
adonis(bc ~ input$map_loaded$Ecosystem)

# PERMDISP
m2 <- betadisper(bc, input$map_loaded$Ecosystem)
anova(m2) # Dispersion not homogeneous



#### Supplementary ####
# Make supplementary table with the metagenome metadata
# Only use the ones used in the analysis - and mark if used for gene, taxa, or both
# Need to get permissions for some of these, also important to mark this (done in Excel)
metadata_gene <- f_meta_MUSiCC %>%
  mutate(Analysis = NA)
metadata_taxa <- input$map_loaded %>%
  select(-count, -GenomeSize) %>%
  mutate(Analysis = NA)

# Names same but formatted differently
names(metadata_gene) <- names(metadata_taxa)

# Combine and fill analysis column with Gene, Taxa, Both
sum(metadata_gene$sampleID %in% metadata_taxa$sampleID) # 2089 both
sum(metadata_taxa$sampleID %notin% metadata_gene$sampleID) # 231 taxa only
sum(metadata_gene$sampleID %notin% metadata_taxa$sampleID) # 22 gene only

d1 <- metadata_gene %>%
  filter(sampleID %in% metadata_taxa$sampleID) %>%
  mutate(Analysis = "Both")

d2 <- metadata_taxa %>%
  filter(sampleID %notin% metadata_gene$sampleID) %>%
  mutate(Analysis = "Taxa")

d3 <- metadata_gene %>%
  filter(sampleID %notin% metadata_taxa$sampleID) %>%
  mutate(Analysis = "Gene")



# lat long got messed up, so reimport it and join it to the used samples, then plot map
lat_long <- read.delim("lyu_habitats_metadata.txt") %>%
  arrange(taxon_oid) %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>%
  separate(., sampleID, into = c("X", "taxon_oid"), sep = 1) %>%
  mutate(taxon_oid = as.numeric(taxon_oid)) %>%
  select(taxon_oid, Latitude, Longitude)

# Asked Christa Pennacchio and she sent me the Gold Database to check restriction status
gold_info <- read_xlsx("~/Desktop/Review/goldData.xlsx", sheet = 6) %>%
  select(`AP IMG TAXON ID`, `AP JGI DATA UTILIZATION STATUS`) %>%
  filter(`AP IMG TAXON ID` %in% metadata_combined$taxon_oid) %>%
  set_names(c("taxon_oid", "JGI Data Utilization Status_GOLD_07Apr2022")) %>%
  mutate(taxon_oid = as.numeric(taxon_oid))
table(gold_info$`JGI Data Utilization Status_GOLD_07Apr2022`)

# Make Table S1
# -Combine and remove unnecessary columns
# -Add accurate lat/long
# Add additional data status column from GOLD
# Reorder columns
metadata_combined <- rbind(d1, d2, d3) %>%
  select(-`IMG Genome ID`, -`Is Public`, -`User Access`, -Salinity, -sampleID, -Latitude, -Longitude) %>%
  left_join(., lat_long, by = "taxon_oid") %>%
  left_join(., gold_info, by = "taxon_oid") %>%
  select(-`Contact Name`, -`Contact Email`, -`JGI Data Utilization Status`, -`JGI Data Utilization Status_GOLD_07Apr2022`, everything(), `Contact Name`, `Contact Email`, `JGI Data Utilization Status`, `JGI Data Utilization Status_GOLD_07Apr2022`)

# Export as Supplementary Table 1
write_xlsx(metadata_combined, "TableS1_Old.xlsx", format_headers = F)

# In Excel make a new column confirming instances in which a PI was contacted and gave permission

# In Excel double check the sample names and remove ones that should be removed
# Save the table and then rerun the analysis filtering by those taxonoids. 



# Make map as Figure S1

# Sample map with ggplot
world <- map_data("world")
ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey", color = "white", size = 0.5) +
  coord_map(projection = "mollweide") +
  labs(x = "Longitude",
       y = "Latitude") + 
  theme_map() +
  geom_point(data = f_meta, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5)

#### ......................... ####
#### ......................... ####

#### REANALYSIS 1 ####
# Rerun analysis after filtering non-relevant sample types (e.g. seawater)
# Also had to remove some restricted status (Kelly Wrighton did not give permission)
# Also only use samples included in both analyses
# Also filter to only metagenomes that have mcrA
# Note: switched hydrogen marker from mtrA to frhA

#### (a) Functional (MUSiCC) ####
# Metadata
f_meta <- read_xlsx("TableS1.xlsx") %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(sampleID = paste("X", taxon_oid, sep = "")) %>%
  arrange(sampleID)

# Functional table (KO)
f <- read.delim("lyu_ko_statanalysis_09-feb-2022/stats_input.txt") %>%
  filter(UniqueID != "GroupID") %>%
  separate(UniqueID, into = c("Text", "KO"), sep = ":") %>%
  select(-Text) %>%
  column_to_rownames(var = "KO") %>%
  replace(is.na(.), 0) %>%
  mutate_if(is.character, as.numeric)

f_comm <- as.data.frame(t(f)) %>%
  filter(rownames(.) %in% f_meta$sampleID) %>%
  arrange(row.names(.)) %>%
  filter(K00399 > 0) # Subset to only metagenomes containing mcrA

# Now export for MUSiCC normalization and correction
for_musicc <- f_comm %>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(var = "KO")
# write.csv(for_musicc, "ko_abundances_updated.csv", row.names = FALSE)
# To run MUSiCC normalization and correction, open Terminal, type jupyter-lab, and run the MUSiCC.ipynb script

# Input
f_comm_MUSiCC <- read.csv("MUSiCC_updated.csv") %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>%
  na.omit()

f_meta_MUSiCC <- f_meta %>%
  filter(sampleID %in% rownames(f_comm_MUSiCC))

# Check IDs
sum(f_meta_MUSiCC$sampleID != rownames(f_comm_MUSiCC))

# Check ecosystem sample size
table(f_meta_MUSiCC$Ecosystem)

# Data frame
gene_plot2_m_lyu <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                               "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                               "MUSiCCmcrA" = f_comm_MUSiCC$K00399,
                               "MUSiCCcdhD" = f_comm_MUSiCC$K00194,
                               "MUSiCCfrhA" = f_comm_MUSiCC$K00440,
                               "MUSiCCmttC" = f_comm_MUSiCC$K14084,
                               "MUSiCCmtbC" = f_comm_MUSiCC$K16179,
                               "MUSiCCmtmC" = f_comm_MUSiCC$K16177,
                               "MUSiCCmtaA" = f_comm_MUSiCC$K14080,
                               "MUSiCCmtsA" = f_comm_MUSiCC$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent"))

#### __ (I) Stats ####
# Run a loop 
kruskal_results_genes <- as.data.frame(matrix(data = NA, 10, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 3:10) {
  k <- kruskal.test(gene_plot2_m_lyu[[i]] ~ gene_plot2_m_lyu$Ecosystem)
  kruskal_results_genes$Gene[i] <- names(gene_plot2_m_lyu)[i]
  kruskal_results_genes$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes$P[i] <- k$p.value
}
kruskal_results_genes <- kruskal_results_genes %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))



#### __ (II) Graph ####
gene_plot2_long_m <- gene_plot2_m_lyu %>%
  pivot_longer(c("MUSiCCmcrA", "MUSiCCcdhD", "MUSiCCfrhA", "MUSiCCmttC", 
                 "MUSiCCmtbC", "MUSiCCmtmC", "MUSiCCmtaA", "MUSiCCmtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "MUSiCCmcrA" = "mcrA",
                              "MUSiCCcdhD" = "cdhD",
                              "MUSiCCfrhA" = "frhA",
                              "MUSiCCmttC" = "mttC",
                              "MUSiCCmtbC" = "mtbC", 
                              "MUSiCCmtmC" = "mtmC", 
                              "MUSiCCmtaA" = "mtaA",
                              "MUSiCCmtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "frhA" = "H2/CO2 (frhA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)"))

lyu_only_m <- gene_plot2_long_m %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

lyu_only_summary_m <- lyu_only_m %>%
  group_by(Hypothesis, Ecosystem, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethylotrophic")

pdf("Figure3_MUSiCCforppt.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (MUSiCC normalized)") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()



#### __ (III) Correlations ####
# Since using only 1 marker gene for each path, check how these are correlated with the other genes in the path
# Check only genes unique to each pathway according to Kurth et al. 2020 Figure 1
check_corrs <- data.frame("mcrA" = f_comm_MUSiCC$K00399,
                          "mcrB" = f_comm_MUSiCC$K00401,
                          "mcrG" = f_comm_MUSiCC$K00402,
                          "cdhD" = f_comm_MUSiCC$K00194,
                          "ackA" = f_comm_MUSiCC$K00925,
                          "pta" = f_comm_MUSiCC$K00625,
                          "eutD" = f_comm_MUSiCC$K04020,
                          "pta2" = f_comm_MUSiCC$K13788,
                          "pta3" = f_comm_MUSiCC$K15024,
                          "acs" = f_comm_MUSiCC$K01895,
                          "cdhC" = f_comm_MUSiCC$K00193,
                          "cdhE" = f_comm_MUSiCC$K00197,
                          "frhA" = f_comm_MUSiCC$K00440,
                          "frhB" = f_comm_MUSiCC$K00441,
                          "frhG" = f_comm_MUSiCC$K00443,
                          "hdrA" = f_comm_MUSiCC$K03388,
                          "hdrB" = f_comm_MUSiCC$K03389,
                          "hdrC" = f_comm_MUSiCC$K03390,
                          # "fdhA" = f_comm_MUSiCC$K22516, not present
                          "fdhB" = f_comm_MUSiCC$K00125,
                          "mvhA" = f_comm_MUSiCC$K14126,
                          "mvhG" = f_comm_MUSiCC$K14128,
                          "mvhD" = f_comm_MUSiCC$K14127,
                          "ftr" = f_comm_MUSiCC$K00672,
                          "mch" = f_comm_MUSiCC$K01499,
                          "hmd" = f_comm_MUSiCC$K13942,
                          "mer" = f_comm_MUSiCC$K00320,
                          "mtrB" = f_comm_MUSiCC$K00578,
                          "mtrC" = f_comm_MUSiCC$K00579,
                          "mtrD" = f_comm_MUSiCC$K00580,
                          "mtrE" = f_comm_MUSiCC$K00581,
                          "mtrF" = f_comm_MUSiCC$K00582,
                          "mtrG" = f_comm_MUSiCC$K00583,
                          "mtrH" = f_comm_MUSiCC$K00584,
                          "mtaA" = f_comm_MUSiCC$K14080,
                          "mtaB" = f_comm_MUSiCC$K04480,
                          "mtaC" = f_comm_MUSiCC$K14081,
                          "mttC" = f_comm_MUSiCC$K14084,
                          "mttB" = f_comm_MUSiCC$K14083,
                          "mtbC" = f_comm_MUSiCC$K16179,
                          "mtbB" = f_comm_MUSiCC$K16178,
                          "mtmC" = f_comm_MUSiCC$K16177,
                          "mtmB" = f_comm_MUSiCC$K16176,
                          "mtbA" = f_comm_MUSiCC$K14082)

# All - mcrA vs mcrB/mcrG
A <- ggplot(check_corrs, aes(mcrA, mcrB)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mcrA", 
       y = "mcrB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

B <- ggplot(check_corrs, aes(mcrA, mcrG)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mcrA", 
       y = "mcrG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Acetate - cdhD vs ackA/pta/eutD/pta2/pta3/acs/cdhC/cdhE
C <- ggplot(check_corrs, aes(cdhD, ackA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "ackA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

D <- ggplot(check_corrs, aes(cdhD, pta)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K00625") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

E <- ggplot(check_corrs, aes(cdhD, eutD)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "eutD") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

f <- ggplot(check_corrs, aes(cdhD, pta2)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K13788") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

G <- ggplot(check_corrs, aes(cdhD, pta3)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K15024") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

H <- ggplot(check_corrs, aes(cdhD, acs)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "acs") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

I <- ggplot(check_corrs, aes(cdhD, cdhC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "cdhC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

J <- ggplot(check_corrs, aes(cdhD, cdhE)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "cdhE") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Hydrogen - frhA vs frhB, frhG, hdrA, hdrB, hdrC, fdhB, mvhA, mvhG, mvhD
K <- ggplot(check_corrs, aes(frhA, frhB)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "frhB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

L <- ggplot(check_corrs, aes(frhA, frhG)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "frhG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

M <- ggplot(check_corrs, aes(frhA, hdrA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

N <- ggplot(check_corrs, aes(frhA, hdrB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

O <- ggplot(check_corrs, aes(frhA, hdrC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

P <- ggplot(check_corrs, aes(frhA, fdhB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "fwdF") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

Q <- ggplot(check_corrs, aes(frhA, fdhB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "fwdG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

R <- ggplot(check_corrs, aes(frhA, mvhA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

S <- ggplot(check_corrs, aes(frhA, mvhG)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

t <- ggplot(check_corrs, aes(frhA, mvhD)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhD") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Methanol - mtaA vs mtaB/mtaC
DD <- ggplot(check_corrs, aes(mtaA, mtaB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtaA", 
       y = "mtaB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

EE <- ggplot(check_corrs, aes(mtaA, mtaC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtaA", 
       y = "mtaC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# TMA - mttC vs mttB/mtbA
FF <- ggplot(check_corrs, aes(mttC, mttB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mttC", 
       y = "mttB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

GG <- ggplot(check_corrs, aes(mttC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mttC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# DMA - mtbC vs mtbB/mtbA
HH <- ggplot(check_corrs, aes(mtbC, mtbB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtbC", 
       y = "mtbB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

II <- ggplot(check_corrs, aes(mtbC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtbC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# MMA - mtmC vs mtmB/mtbA
JJ <- ggplot(check_corrs, aes(mtmC, mtmB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtmC", 
       y = "mtmB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

KK <- ggplot(check_corrs, aes(mtmC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtmC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# DMS/MeSH/MMPA - mtsA vs mtsB (cannot do because no KO)

# Plot (adjust device size and save manually)
plot_grid(A, B, NULL, NULL, NULL,
          C, D, E, f, G, 
          H, I, J, NULL, NULL,
          K, L, M, N, NULL,
          O, P, Q, R, S,
          DD, EE, NULL, FF, GG,
          HH, II, NULL, JJ, KK,
          ncol = 5, nrow = 7, 
          align = "hv")



#### (b) Taxonomic (family) ####
# Import
tax_table_fp <- file.path("~/Desktop/Review/IMG_families_mctoolsr.txt")
map_fp <- file.path("~/Desktop/Review/lyu_habitats_metadata_updated.txt")
input = load_taxa_table(tax_table_fp, map_fp)
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""))

# Match the functional analysis
input = filter_data(input,
                    filter_cat = "sampleID",
                    keep_vals = f_meta_MUSiCC$sampleID)

# Check sequencing depth 
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded))
se(colSums(input$data_loaded))
input$map_loaded$count <- colSums(input$data_loaded)
ggplot(input$map_loaded, aes(reorder(Ecosystem, count, mean), count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "# Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

ggplot(input$map_loaded, aes(`Genome Size   * assembled`, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  labs(x = "Assembled genome size", 
       y = "Assigned family reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# Rename size column
input$map_loaded$GenomeSize = input$map_loaded$`Genome Size   * assembled`

# Phyla
tax_sum_phyla <- summarize_taxonomy(input, level = 2, report_higher_tax = T, relative = FALSE)
# Note there are no Bathyarchaeota or Verstraetearchaeota

# Summarize by family, extract methanogens, calculate CPM
tax_sum_family_wTax <- summarize_taxonomy(input, level = 5, report_higher_tax = T, relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
methano_wTax <- methano_wTax[!grepl("Plasmid", rownames(methano_wTax)),]
# Note this also contains a methanotroph Methermicoccaceae and a plasmid Methanobacteriales

# Without higher tax and no plasmid
tax_sum_family <- summarize_taxonomy(input, level = 5, report_higher_tax = FALSE, relative = FALSE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
methano <- subset(methano, rownames(methano) != "Methanobacteriales")
tax_sum_family_df <- t(methano) %>%
  as.data.frame() %>%
  mutate(`Candidatus Methanoperedenaceae` = (`Candidatus Methanoperedenaceae`*1000000)/input$map_loaded$GenomeSize,
         Methanobacteriaceae = (Methanobacteriaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocaldococcaceae = (Methanocaldococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocellaceae = (Methanocellaceae*1000000)/input$map_loaded$GenomeSize,
         Methanococcaceae = (Methanococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocorpusculaceae = (Methanocorpusculaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomassiliicoccaceae = (Methanomassiliicoccaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomicrobiaceae = (Methanomicrobiaceae*1000000)/input$map_loaded$GenomeSize,
         Methanonatronarchaeaceae = (Methanonatronarchaeaceae*1000000)/input$map_loaded$GenomeSize,
         Methanopyraceae = (Methanopyraceae*1000000)/input$map_loaded$GenomeSize,
         Methanoregulaceae = (Methanoregulaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosaetaceae = (Methanosaetaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosarcinaceae = (Methanosarcinaceae*1000000)/input$map_loaded$GenomeSize,
         Methanospirillaceae = (Methanospirillaceae*1000000)/input$map_loaded$GenomeSize,
         Methanothermaceae = (Methanothermaceae*1000000)/input$map_loaded$GenomeSize,
         Methanotrichaceae = (Methanotrichaceae*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanomicrobiales` = (`unclassified Methanomicrobiales`*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanosarcinales` = (`unclassified Methanosarcinales`*1000000)/input$map_loaded$GenomeSize) 

input$map_loaded <- input$map_loaded %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))
sum(input$map$sampleID != rownames(tax_sum_family_df))
input$map_loaded <- cbind(input$map_loaded, tax_sum_family_df)



#### __(I) Stats ####
# Run a loop
kruskal_results_taxa <- as.data.frame(matrix(data = NA, 47, 3)) %>%
  set_names(c("Family", "X2", "P"))
for (i in 30:47) {
  k <- kruskal.test(input$map_loaded[[i]] ~ input$map_loaded$Ecosystem)
  kruskal_results_taxa$Family[i] <- names(input$map_loaded)[i]
  kruskal_results_taxa$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_taxa$P[i] <- k$p.value
}
kruskal_results_taxa <- kruskal_results_taxa %>%
  filter(is.na(Family) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))


#### __(II) Graph ####
family_plot_long <- pivot_longer(input$map_loaded,
                                 cols = names(input$map_loaded)[30:47],
                                 names_to = "Family",
                                 values_to = "Abundance")
test <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  group_by(Family) %>%
  summarise(max_abund = max(mean_abund))

bars <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  filter(Family != "unclassified Methanomicrobiales",
         Family != "unclassified Methanosarcinales",
         Family != "Candidatus Methanoperedenaceae")

filter_bars <- bars %>%
  group_by(Family) %>%
  summarise(max_abund = max(mean_abund)) %>%
  filter(max_abund > 1)

bars_top <- bars %>%
  filter(Family %in% filter_bars$Family)

sort_bars <- bars_top %>%
  group_by(Family) %>%
  summarise(sum_abund = sum(mean_abund)) %>%
  arrange(sum_abund)

bars_top_sorted <- bars_top %>%
  mutate(Family = factor(Family,
                         levels = sort_bars$Family))

bars_top_sorted_guild <- bars_top_sorted %>%
  mutate(Family = factor(Family,
                         levels = c("Methanosaetaceae", "Methanotrichaceae",
                                    "Methanospirillaceae", "Methanocellaceae",
                                    "Methanocorpusculaceae", "Methanococcaceae",
                                    "Methanomicrobiaceae", "Methanocaldococcaceae",
                                    "Methanobacteriaceae", "Methanoregulaceae",
                                    "Methanomassiliicoccaceae", "Methanosarcinaceae")))

n_families <- bars %>%
  group_by(Hypothesis, Ecosystem) %>%
  summarise(n_fam = sum(mean_abund > 0),
            tot = sum(mean_abund)) %>%
  mutate(y = tot + 2)

nb.cols <- nrow(filter_bars)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethylotrophic")

# Figure 4 - finalize in powerpoint
pdf("Figure4_forppt.pdf", width = 8.5, height = 4.5)
ggplot(bars_top_sorted_guild, aes(Ecosystem, mean_abund, fill = Family)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.box.margin=margin(0,-5,0,15))
dev.off()

# Remake figure aggregated by guild
bars_top_sorted_guild$Guild <- NA
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanosaetaceae") {
    bars_top_sorted_guild$Guild[i] <- "Acetoclasts"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanotrichaceae") {
    bars_top_sorted_guild$Guild[i] <- "Acetoclasts"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanospirillaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocellaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocorpusculaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanococcaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanomicrobiaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocaldococcaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanobacteriaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanoregulaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophs"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanomassiliicoccaceae") {
    bars_top_sorted_guild$Guild[i] <- "Methylotrophs (H2-dep.)"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanosarcinaceae") {
    bars_top_sorted_guild$Guild[i] <- "Mixotrophs"
  }
}

bars_top_sorted_guild_ag <- bars_top_sorted_guild %>%
  group_by(Ecosystem, Hypothesis, Guild) %>%
  summarise(mean_abund = sum(mean_abund))
  

# Guild figure
pdf("FigureS3.pdf", width = 8.5, height = 4.5)
ggplot(bars_top_sorted_guild_ag, aes(Ecosystem, mean_abund, fill = Guild)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  scale_fill_discrete(breaks = bars_top_sorted_guild_ag$Guild,
                      labels = c("Acetoclasts",
                                 "Hydrogenotrophs",
                                 bquote('Methylotrophs ('*H[2]~'dep.)'),
                                 "Mixotrophs")) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7))
dev.off()

png("FigureS3.png", width = 8.5, height = 4.5, unit = "in", res = 300)
ggplot(bars_top_sorted_guild_ag, aes(Ecosystem, mean_abund, fill = Guild)) +
  geom_bar(stat = "identity", colour = "black", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  scale_fill_discrete(breaks = bars_top_sorted_guild_ag$Guild,
                      labels = c("Acetoclasts",
                                 "Hydrogenotrophs",
                                 bquote('Methylotrophs ('*H[2]~'dep.)'),
                                 "Mixotrophs")) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7))
dev.off()

hue_pal()(4)

# Sample size
table(input$map_loaded$Ecosystem)



#### Map ####
# Sample map with ggplot
world <- map_data("world")

pdf("FigureS1.pdf", width = 7, height = 5)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  geom_point(data = f_meta, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

png("FigureS1.png", width = 7, height = 5, units = "in", res = 300)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  geom_point(data = f_meta, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



#### Table S1 ####
# Save Updated Table S1 with the 395 samples actually analyzed
# -Combine and remove unnecessary columns
# -Add additional data status column from GOLD
# -Reorder columns
metadata <- f_meta_MUSiCC %>%
  select(-sampleID) %>%
  select(-`Contact Name`, -`Contact Email`, -`JGI Data Utilization Status`, -`JGI Data Utilization Status_GOLD_07Apr2022`, -PermissionNotes, everything(), `Contact Name`, `Contact Email`, `JGI Data Utilization Status`, `JGI Data Utilization Status_GOLD_07Apr2022`, PermissionNotes)
names(metadata)

# Export as Supplementary Table 1
# write_xlsx(metadata, "TableS1_Updated_395.xlsx", format_headers = F)

# In Excel fill PermissionNotes column confirming instances in which a PI was contacted and gave permission

#### .......................... ####
#### .......................... ####
#### REANALYSIS 2 ####
#### Add coastal/hypersaline ####
# Get coastal and hypersaline samples, with mcrA, with permission
# Then add to f_comm and f_meta, re-normalize with MUSiCC, reanalyze
tableS1 <- read_xlsx("~/Desktop/Review/TableS1_Updated_395.xlsx")

f_meta <- readRDS("f_meta.rds") %>%
  mutate(sampleID = paste("X", IMG.Genome.ID, sep = ""))

f_meta_hypcoast <- subset(f_meta, Ecosystem == "Hypersaline" | 
                            Ecosystem == "Wetland (coastal)")

gold_info <- read_xlsx("~/Desktop/Review/goldData.xlsx", sheet = 6) %>%
  select(`AP IMG TAXON ID`, `AP JGI DATA UTILIZATION STATUS`) %>%
  filter(`AP IMG TAXON ID` %in% f_meta_hypcoast$taxon_oid) %>%
  set_names(c("taxon_oid", "JGI Data Utilization Status_GOLD_07Apr2022")) %>%
  mutate(taxon_oid = as.numeric(taxon_oid))

f_comm <- readRDS("f_comm.rds")
f_meta <- readRDS("f_meta.rds") %>%
  mutate(sampleID = paste("X", IMG.Genome.ID, sep = ""))

f_meta_hypcoast <- subset(f_meta, Ecosystem == "Hypersaline" | 
                            Ecosystem == "Wetland (coastal)") %>%
  left_join(., gold_info, by = "taxon_oid")

f_comm_hypcoast <- subset(f_comm, rownames(f_comm) %in% f_meta_hypcoast$sampleID)

# With mcrA
f_comm_hypcoast <- subset(f_comm_hypcoast, K00399 > 0)
f_meta_hypcoast <- subset(f_meta_hypcoast, f_meta_hypcoast$sampleID %in% row.names(f_comm_hypcoast))

# Not already in Table S1_Updated_395
f_meta_hypcoast <- subset(f_meta_hypcoast, f_meta_hypcoast$taxon_oid %notin% tableS1$taxon_oid)
f_comm_hypcoast <- subset(f_comm, rownames(f_comm) %in% f_meta_hypcoast$sampleID)

# With permission (unrestricted or our projects)
f_meta_hypcoast <- subset(f_meta_hypcoast, `JGI Data Utilization Status_GOLD_07Apr2022` != "Restricted" |
                            Study.Name == "Natural and restored wetland microbial communities from the San Francisco Bay, California, USA, that impact long-term carbon sequestration")
f_comm_hypcoast <- subset(f_comm_hypcoast, rownames(f_comm_hypcoast) %in% f_meta_hypcoast$sampleID)

# Now append to the previous 395 and rerun
f_meta_hypcoast <- f_meta_hypcoast %>%
  mutate(PermissionNotes = "NA") %>%
  select(taxon_oid, Ecosystem, -Domain, Study.Name, Genome.Name...Sample.Name,
         Sequencing.Center, -IMG.Genome.ID, Add.Date, Assembly.Method, Has.Coverage,
         -Is.Public, Release.Date, -User.Access, Sequencing.Method,
         Genome.Size.....assembled, Gene.Count.....assembled, GC.....assembled,
         KO.......assembled, richness_KO, Hypothesis, Latitude, Longitude,
         Contact.Email, Contact.Name, JGI.Data.Utilization.Status,
         `JGI Data Utilization Status_GOLD_07Apr2022`, PermissionNotes, sampleID)

#### (a) Functional (MUSiCC) ####
# Metadata
f_meta <- read_xlsx("TableS1.xlsx") %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(sampleID = paste("X", taxon_oid, sep = "")) %>%
  arrange(sampleID)

f_meta_hypcoast <- f_meta_hypcoast %>%
  set_names(names(f_meta))

f_meta <- rbind(f_meta, f_meta_hypcoast)

# Functional table (KO)
f <- read.delim("lyu_ko_statanalysis_09-feb-2022/stats_input.txt") %>%
  filter(UniqueID != "GroupID") %>%
  separate(UniqueID, into = c("Text", "KO"), sep = ":") %>%
  select(-Text) %>%
  column_to_rownames(var = "KO") %>%
  replace(is.na(.), 0) %>%
  mutate_if(is.character, as.numeric)

f_comm <- as.data.frame(t(f)) %>%
  filter(rownames(.) %in% f_meta$sampleID) %>%
  arrange(row.names(.)) %>%
  filter(K00399 > 0) # Subset to only metagenomes containing mcrA

# Now export for MUSiCC normalization and correction
for_musicc <- f_comm %>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric) %>%
  rownames_to_column(var = "KO")
# write.csv(for_musicc, "ko_abundances_updated2.csv", row.names = FALSE)
# To run MUSiCC normalization and correction, open Terminal, type jupyter-lab, and run the MUSiCC.ipynb script

# Input
f_comm_MUSiCC <- read.csv("MUSiCC_updated2.csv") %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>%
  na.omit() %>%
  arrange(rownames(.))

f_meta_MUSiCC <- f_meta %>%
  filter(sampleID %in% rownames(f_comm_MUSiCC)) %>%
  arrange(sampleID)

# n = 465. Some dropped because NA in MUSiCC
# Previously input n = 397, output n = 395
# Now input n = 489 (because +92 coastal/hypersaline added), n output = 465
# So total added is 70

# Check IDs
sum(f_meta_MUSiCC$sampleID != rownames(f_comm_MUSiCC))

# Check ecosystem sample size
table(f_meta_MUSiCC$Ecosystem)

# Save!
# saveRDS(f_meta_MUSiCC, "f_meta_MUSiCC_reanalysis2_n465.rds")
# saveRDS(f_comm_MUSiCC, "f_comm_MUSiCCreanalysis2_n465.rds")
f_meta_MUSiCC <- readRDS("f_meta_MUSiCC_reanalysis2_n465.rds")
f_comm_MUSiCC <- readRDS("f_comm_MUSiCCreanalysis2_n465.rds")


# Data frame
gene_plot2_m_lyu <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                               "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                               "MUSiCCmcrA" = f_comm_MUSiCC$K00399,
                               "MUSiCCcdhD" = f_comm_MUSiCC$K00194,
                               "MUSiCCfrhA" = f_comm_MUSiCC$K00440,
                               "MUSiCCmttC" = f_comm_MUSiCC$K14084,
                               "MUSiCCmtbC" = f_comm_MUSiCC$K16179,
                               "MUSiCCmtmC" = f_comm_MUSiCC$K16177,
                               "MUSiCCmtaA" = f_comm_MUSiCC$K14080,
                               "MUSiCCmtsA" = f_comm_MUSiCC$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))
levels(gene_plot2_m_lyu$Ecosystem)

#### __ (I) Stats ####
# Run a loop 
kruskal_results_genes <- as.data.frame(matrix(data = NA, 10, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 3:10) {
  k <- kruskal.test(gene_plot2_m_lyu[[i]] ~ gene_plot2_m_lyu$Ecosystem)
  kruskal_results_genes$Gene[i] <- names(gene_plot2_m_lyu)[i]
  kruskal_results_genes$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes$P[i] <- k$p.value
}
kruskal_results_genes <- kruskal_results_genes %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))



#### __ (II) Graph ####
gene_plot2_long_m <- gene_plot2_m_lyu %>%
  pivot_longer(c("MUSiCCmcrA", "MUSiCCcdhD", "MUSiCCfrhA", "MUSiCCmttC", 
                 "MUSiCCmtbC", "MUSiCCmtmC", "MUSiCCmtaA", "MUSiCCmtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "MUSiCCmcrA" = "mcrA",
                              "MUSiCCcdhD" = "cdhD",
                              "MUSiCCfrhA" = "frhA",
                              "MUSiCCmttC" = "mttC",
                              "MUSiCCmtbC" = "mtbC", 
                              "MUSiCCmtmC" = "mtmC", 
                              "MUSiCCmtaA" = "mtaA",
                              "MUSiCCmtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "frhA" = "H2/CO2 (frhA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)"))

lyu_only_m <- gene_plot2_long_m %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

lyu_only_summary_m <- lyu_only_m %>%
  group_by(Hypothesis, Ecosystem, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethyl-based")

pdf("Figure4_MUSiCCforppt.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (MUSiCC normalized)") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()



#### __ (III) Correlations ####
# Since using only 1 marker gene for each path, check how these are correlated with the other genes in the path
# Check only genes unique to each pathway according to Kurth et al. 2020 Figure 1
check_corrs <- data.frame("mcrA" = f_comm_MUSiCC$K00399,
                          "mcrB" = f_comm_MUSiCC$K00401,
                          "mcrG" = f_comm_MUSiCC$K00402,
                          "cdhD" = f_comm_MUSiCC$K00194,
                          "ackA" = f_comm_MUSiCC$K00925,
                          "pta" = f_comm_MUSiCC$K00625,
                          "eutD" = f_comm_MUSiCC$K04020,
                          "pta2" = f_comm_MUSiCC$K13788,
                          "pta3" = f_comm_MUSiCC$K15024,
                          "acs" = f_comm_MUSiCC$K01895,
                          "cdhC" = f_comm_MUSiCC$K00193,
                          "cdhE" = f_comm_MUSiCC$K00197,
                          "frhA" = f_comm_MUSiCC$K00440,
                          "frhB" = f_comm_MUSiCC$K00441,
                          "frhG" = f_comm_MUSiCC$K00443,
                          "hdrA" = f_comm_MUSiCC$K03388,
                          "hdrB" = f_comm_MUSiCC$K03389,
                          "hdrC" = f_comm_MUSiCC$K03390,
                          # "fdhA" = f_comm_MUSiCC$K22516, not present
                          "fdhB" = f_comm_MUSiCC$K00125,
                          "mvhA" = f_comm_MUSiCC$K14126,
                          "mvhG" = f_comm_MUSiCC$K14128,
                          "mvhD" = f_comm_MUSiCC$K14127,
                          "ftr" = f_comm_MUSiCC$K00672,
                          "mch" = f_comm_MUSiCC$K01499,
                          "hmd" = f_comm_MUSiCC$K13942,
                          "mer" = f_comm_MUSiCC$K00320,
                          "mtrB" = f_comm_MUSiCC$K00578,
                          "mtrC" = f_comm_MUSiCC$K00579,
                          "mtrD" = f_comm_MUSiCC$K00580,
                          "mtrE" = f_comm_MUSiCC$K00581,
                          "mtrF" = f_comm_MUSiCC$K00582,
                          "mtrG" = f_comm_MUSiCC$K00583,
                          "mtrH" = f_comm_MUSiCC$K00584,
                          "mtaA" = f_comm_MUSiCC$K14080,
                          "mtaB" = f_comm_MUSiCC$K04480,
                          "mtaC" = f_comm_MUSiCC$K14081,
                          "mttC" = f_comm_MUSiCC$K14084,
                          "mttB" = f_comm_MUSiCC$K14083,
                          "mtbC" = f_comm_MUSiCC$K16179,
                          "mtbB" = f_comm_MUSiCC$K16178,
                          "mtmC" = f_comm_MUSiCC$K16177,
                          "mtmB" = f_comm_MUSiCC$K16176,
                          "mtbA" = f_comm_MUSiCC$K14082)

# All - mcrA vs mcrB/mcrG
A <- ggplot(check_corrs, aes(mcrA, mcrB)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mcrA", 
       y = "mcrB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

B <- ggplot(check_corrs, aes(mcrA, mcrG)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mcrA", 
       y = "mcrG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Acetate - cdhD vs ackA/pta/eutD/pta2/pta3/acs/cdhC/cdhE
C <- ggplot(check_corrs, aes(cdhD, ackA)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "ackA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

D <- ggplot(check_corrs, aes(cdhD, pta)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K00625") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

E <- ggplot(check_corrs, aes(cdhD, eutD)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "eutD") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

f <- ggplot(check_corrs, aes(cdhD, pta2)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K13788") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

G <- ggplot(check_corrs, aes(cdhD, pta3)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "pta K15024") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

H <- ggplot(check_corrs, aes(cdhD, acs)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "acs") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

I <- ggplot(check_corrs, aes(cdhD, cdhC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "cdhC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

J <- ggplot(check_corrs, aes(cdhD, cdhE)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "cdhD", 
       y = "cdhE") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Hydrogen - frhA vs frhB, frhG, hdrA, hdrB, hdrC, fdhB, mvhA, mvhG, mvhD
K <- ggplot(check_corrs, aes(frhA, frhB)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "frhB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

L <- ggplot(check_corrs, aes(frhA, frhG)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "frhG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

M <- ggplot(check_corrs, aes(frhA, hdrA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

N <- ggplot(check_corrs, aes(frhA, hdrB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

O <- ggplot(check_corrs, aes(frhA, hdrC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "hdrC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

P <- ggplot(check_corrs, aes(frhA, fdhB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "fwdF") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

Q <- ggplot(check_corrs, aes(frhA, fdhB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "fwdG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

R <- ggplot(check_corrs, aes(frhA, mvhA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

S <- ggplot(check_corrs, aes(frhA, mvhG)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_point(size = 1, alpha = 0.1) +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhG") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

t <- ggplot(check_corrs, aes(frhA, mvhD)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "frhA", 
       y = "mvhD") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# Methanol - mtaA vs mtaB/mtaC
DD <- ggplot(check_corrs, aes(mtaA, mtaB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtaA", 
       y = "mtaB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

EE <- ggplot(check_corrs, aes(mtaA, mtaC)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtaA", 
       y = "mtaC") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# TMA - mttC vs mttB/mtbA
FF <- ggplot(check_corrs, aes(mttC, mttB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mttC", 
       y = "mttB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

GG <- ggplot(check_corrs, aes(mttC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mttC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# DMA - mtbC vs mtbB/mtbA
HH <- ggplot(check_corrs, aes(mtbC, mtbB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtbC", 
       y = "mtbB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

II <- ggplot(check_corrs, aes(mtbC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtbC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# MMA - mtmC vs mtmB/mtbA
JJ <- ggplot(check_corrs, aes(mtmC, mtmB)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtmC", 
       y = "mtmB") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

KK <- ggplot(check_corrs, aes(mtmC, mtbA)) +
  geom_point(size = 1, alpha = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
  geom_smooth(method = "lm", size = 0.2) +
  labs(x = "mtmC", 
       y = "mtbA") +
  theme_bw() +
  theme(axis.title = element_text(size = 8, face = "italic"),
        axis.text = element_text(size = 6))

# DMS/MeSH/MMPA - mtsA vs mtsB (cannot do because no KO)

# Plot (adjust device size and save manually)
plot_grid(A, B, NULL, NULL, NULL,
          C, D, E, f, G, 
          H, I, J, NULL, NULL,
          K, L, M, N, NULL,
          O, P, Q, R, S,
          DD, EE, NULL, FF, GG,
          HH, II, NULL, JJ, KK,
          ncol = 5, nrow = 7, 
          align = "hv")



#### (aa) Functional (Split) ####
# Reanalyze KO abundances but for just the archaeal proportion of the metagenomes
# Sent Dongying Wu of IMG staff list of taxonoids and taxa "Archaea" and "Bacteria"
# Dongying ran custom python scripts on JGI super computer to pull out KOs of only Archaea and Bacteria
# Folder Archaea_KO_465 has a file for each metagenome with the KO hits of the Archaea scaffolds
# Folder Bacteria_KO_465 has a file for each metagenome with the KO hits of the Bacteria scaffolds
# Then subtract these from the full KO profile to get "Other"

#### __ (I) DESeq ####
# Run a for loop to read in the file for each metagenome and combine into 1 (takes a while)
setwd("~/Desktop/Review/Archaea_KO_465/")
list.files()
ko <- list()
ko_input <- data.frame(V1 = "NA",
                       V2 = "NA",
                       V3 = "NA")
ko_table <- ko_input
for (i in 1:length(list.files())) {
  ko[[i]] <- read.delim(list.files()[i], header = F)
  ko_table <- ko_table %>%
    rbind(ko[[i]])
}
setwd("~/Desktop/Review")

# Clean up table (takes a while, ~14 minutes)
ko_table_wTax <- ko_table %>%
  filter(V1 != "NA") %>%
  separate(V1, into = c("Junk", "KO"), sep = ":") %>%
  dplyr::select(-Junk, -V2) %>%
  mutate(taxon_oid = substring(V3, first = 1, last = 10)) %>%
  separate(V3, into = c("Junk", "taxonomy"), sep = "Archaea;") %>%
  dplyr::select(-Junk) %>%
  separate(taxonomy, 
           into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>%
  dplyr::select(taxon_oid, KO, everything())

# KO count (abundance) by metagenome
ko_table_MGcount <- ko_table_wTax %>%
  dplyr::select(taxon_oid, KO) %>%
  group_by(taxon_oid, KO) %>%
  summarise(KO_abund = n()) %>%
  pivot_wider(id_cols = KO, names_from = taxon_oid, values_from = KO_abund) %>%
  column_to_rownames(var = "KO")
ko_table_MGcount[is.na(ko_table_MGcount) == TRUE] <- 0

# Make community style table and metadata, match IDs
ko_comm <- ko_table_MGcount %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% f_meta_MUSiCC$taxon_oid) %>%
  arrange(rownames(.))
sum(rownames(ko_comm) != f_meta_MUSiCC$taxon_oid)
rownames(ko_comm) <- f_meta_MUSiCC$sampleID
sum(rownames(ko_comm) != f_meta_MUSiCC$sampleID)
rs <- as.data.frame(rowSums(ko_comm))
sort(-colSums(ko_comm))

# Now export for MUSiCC normalization and correction
#for_musicc <- ko_comm %>%
#  t() %>%
#  as.data.frame() %>%
#  mutate_if(is.character, as.numeric) %>%
#  rownames_to_column(var = "KO")
# write.csv(for_musicc, "ko_abundances_archaea.csv", row.names = FALSE)
# To run MUSiCC normalization and correction, open Terminal, type jupyter-lab, and run the MUSiCC.ipynb script
# Getting errors, test a subset
#test <- for_musicc[,c(1:57, 59:179, 181:466)]
#write.csv(test, "test.csv", row.names = FALSE)
# Also many without USiCGs
# Perhaps DESeq better for this? Or CPM? Of archaeal or assembled reads?

# Input
#f_comm_MUSiCC_a <- read.csv("MUSiCC_archaea.csv") %>%
#  column_to_rownames(var = "KO") %>%
#  t() %>%
#  as.data.frame() %>%
#  na.omit() %>%
#  arrange(rownames(.))

#ko_meta <- input_fungi$map_loaded %>%
#  filter(taxon_oid %in% rownames(ko_comm)) %>%
#  arrange(taxon_oid) %>%
#  left_join(., methods, by = "sampleID") %>%
#  mutate(rn = sampleID) %>%
#  column_to_rownames(var = "rn")

# DESeq Normalization (run once, then reload from saved)
#dds <- DESeqDataSetFromMatrix(countData = t(ko_comm) + 1,
#                              colData = f_meta_MUSiCC,
#                              design = ~ 1)
#dds <- estimateSizeFactors(dds)
#dds <- estimateDispersions(dds)
#ko_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))
# Save so you don't have to redo the DESeq (takes a while)
# saveRDS(ko_comm_DESeq, "ko_comm_DESeq.rds")
ko_comm_DESeq <- readRDS("ko_comm_DESeq.rds")

# Data frame
gene_plot2_m_lyu <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                               "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                               "MUSiCCmcrA" = ko_comm_DESeq$K00399,
                               "MUSiCCcdhD" = ko_comm_DESeq$K00194,
                               "MUSiCCfrhA" = ko_comm_DESeq$K00440,
                               "MUSiCCmttC" = ko_comm_DESeq$K14084,
                               "MUSiCCmtbC" = ko_comm_DESeq$K16179,
                               "MUSiCCmtmC" = ko_comm_DESeq$K16177,
                               "MUSiCCmtaA" = ko_comm_DESeq$K14080,
                               "MUSiCCmtsA" = ko_comm_DESeq$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))
levels(gene_plot2_m_lyu$Ecosystem)

gene_plot2_long_m <- gene_plot2_m_lyu %>%
  pivot_longer(c("MUSiCCmcrA", "MUSiCCcdhD", "MUSiCCfrhA", "MUSiCCmttC", 
                 "MUSiCCmtbC", "MUSiCCmtmC", "MUSiCCmtaA", "MUSiCCmtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "MUSiCCmcrA" = "mcrA",
                              "MUSiCCcdhD" = "cdhD",
                              "MUSiCCfrhA" = "frhA",
                              "MUSiCCmttC" = "mttC",
                              "MUSiCCmtbC" = "mtbC", 
                              "MUSiCCmtmC" = "mtmC", 
                              "MUSiCCmtaA" = "mtaA",
                              "MUSiCCmtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "frhA" = "H2/CO2 (frhA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)"))

lyu_only_m <- gene_plot2_long_m %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

lyu_only_summary_m <- lyu_only_m %>%
  group_by(Hypothesis, Ecosystem, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethyl-based")

pdf("Figure4_DESeqforppt_arch.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "Abundance (DESeq normalized)") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()


#### __ (II) CPM ####
# Get KO CPM for Archael, Bacterial, Other

# KOs
ko_list <- c("KO:K00399", "KO:K00194", "KO:K00440", "KO:K14084", 
             "KO:K16179", "KO:K16177", "KO:K14080", "KO:K16954")
ko_list2 <- c("K00399", "K00194", "K00440", "K14084", 
             "K16179", "K16177", "K14080", "K16954")

### Archaea
# Run a for loop to read in the file for each metagenome and combine into 1 (takes a while)
setwd("~/Desktop/Review/Archaea_KO_465/")
list.files()
ko_a <- list()
ko_input_a <- data.frame(V1 = "NA",
                         V2 = "NA",
                         V3 = "NA")
ko_table_a <- ko_input_a
for (i in 1:length(list.files())) {
  ko_a[[i]] <- read.delim(list.files()[i], header = F)
  ko_table_a <- ko_table_a %>%
    rbind(ko_a[[i]])
}

# Clean up table (takes a while, ~14 minutes)
ko_table_wTax_a <- ko_table_a %>%
  filter(V1 %in% ko_list) %>%
  separate(V1, into = c("Junk", "KO"), sep = ":") %>%
  dplyr::select(-Junk, -V2) %>%
  mutate(taxon_oid = substring(V3, first = 1, last = 10)) %>%
  separate(V3, into = c("Junk", "taxonomy"), sep = "Archaea;") %>%
  dplyr::select(-Junk) %>%
  separate(taxonomy, 
           into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>%
  dplyr::select(taxon_oid, KO, everything())

# KO count (abundance) by metagenome
ko_table_MGcount_a <- ko_table_wTax_a %>%
  dplyr::select(taxon_oid, KO) %>%
  group_by(taxon_oid, KO) %>%
  summarise(KO_abund = n()) %>%
  pivot_wider(id_cols = KO, names_from = taxon_oid, values_from = KO_abund) %>%
  column_to_rownames(var = "KO")
ko_table_MGcount_a[is.na(ko_table_MGcount_a) == TRUE] <- 0

# Add any metagenomes that are all zeroes
missing_a <- f_meta_MUSiCC %>%
  filter(f_meta_MUSiCC$taxon_oid %notin% colnames(ko_table_MGcount_a)) %>%
  dplyr::select(taxon_oid)

ko_table_MGcount_a <- ko_table_MGcount_a %>%
  mutate("3300002586" = 0,
         "3300027887" = 0)

# Make community style table and metadata, match IDs
ko_comm_a <- ko_table_MGcount_a %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% f_meta_MUSiCC$taxon_oid) %>%
  arrange(rownames(.))
sum(rownames(ko_comm_a) != f_meta_MUSiCC$taxon_oid)
# saveRDS(ko_comm_a, "~/Desktop/Review/ko_comm_a.rds")
ko_comm_a <- readRDS("ko_comm_a.rds")



### Bacteria
# Run a for loop to read in the file for each metagenome and combine into 1 (takes a while)
# PC crashing, files too big, Do this part on Cori supercomputer!!
# PrepBacterialKOTable.ipynb
setwd("~/Desktop/Review/Bacteria_KO_465/")
list.files()
ko_b <- list()
ko_input_b <- data.frame(V1 = "NA",
                         V2 = "NA",
                         V3 = "NA")
ko_table_b <- ko_input_b
for (i in 93:100) {
  ko_b[[i]] <- read.delim(list.files()[i], header = F)
  ko_table_b <- ko_table_b %>%
    rbind(ko_b[[i]])
}

# Clean up table
ko_table_wTax_b <- ko_table_b %>%
  filter(V1 %in% ko_list) %>%
  filter(V1 != "NA") %>%
  separate(V1, into = c("Junk", "KO"), sep = ":") %>%
  dplyr::select(-Junk, -V2) %>%
  mutate(taxon_oid = substring(V3, first = 1, last = 10)) %>%
  separate(V3, into = c("Junk", "taxonomy"), sep = "Archaea;") %>%
  dplyr::select(-Junk) %>%
  separate(taxonomy, 
           into = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = ";") %>%
  dplyr::select(taxon_oid, KO, everything())

# KO count (abundance) by metagenome
ko_table_MGcount_b <- ko_table_wTax_b %>%
  dplyr::select(taxon_oid, KO) %>%
  group_by(taxon_oid, KO) %>%
  summarise(KO_bbund = n()) %>%
  pivot_wider(id_cols = KO, names_from = taxon_oid, values_from = KO_bbund) %>%
  column_to_rownames(var = "KO")
ko_table_MGcount_b[is.na(ko_table_MGcount_b) == TRUE] <- 0

# Add any metagenomes that are all zeroes
missing_b <- f_meta_MUSiCC %>%
  filter(f_meta_MUSiCC$taxon_oid %notin% colnames(ko_table_MGcount_b)) %>%
  dplyr::select(taxon_oid)

ko_table_MGcount_b <- ko_table_MGcount_b %>%
  mutate("3300002586" = 0,
         "3300027887" = 0)

# Make community style table and metadata, match IDs
ko_comm_b <- ko_table_MGcount_b %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% f_meta_MUSiCC$taxon_oid) %>%
  arrange(rownames(.))
sum(rownames(ko_comm_b) != f_meta_MUSiCC$taxon_oid)
ko_comm_b <- readRDS("ko_comm_b.rds") %>%
  dplyr::select(colnames(ko_comm_a))

### Other
# Get full table and subtract bacteria and archaea
f <- read.delim("lyu_ko_statanalysis_09-feb-2022/stats_input.txt") %>%
  filter(UniqueID != "GroupID") %>%
  separate(UniqueID, into = c("Text", "KO"), sep = ":") %>%
  select(-Text) %>%
  column_to_rownames(var = "KO") %>%
  replace(is.na(.), 0) %>%
  mutate_if(is.character, as.numeric)
ko_comm_t <- f %>%
  filter(rownames(.) %in% ko_list2) %>%
  set_names(substr(colnames(.), 2, 11)) %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% f_meta_MUSiCC$taxon_oid) %>%
  arrange(rownames(.))
# saveRDS(ko_comm_t, "ko_comm_t.rds")
ko_comm_t <- readRDS("ko_comm_t.rds")

# Matrix subtraction
# Make sure rownames and column names match!
sum(colnames(ko_comm_a) != colnames(ko_comm_b))
sum(colnames(ko_comm_a) != colnames(ko_comm_t))
mat_a <- as.matrix(ko_comm_a)
mat_b <- as.matrix(ko_comm_b)
mat_t <- as.matrix(ko_comm_t)
mat_o <- mat_t - mat_a - mat_b
ko_comm_o <- as.data.frame(mat_o)
# saveRDS(ko_comm_o, "ko_comm_o.rds")
ko_comm_o <- readRDS("ko_comm_o.rds")

#### ____Load ####
f_meta_MUSiCC <- readRDS("f_meta_MUSiCC_reanalysis2_n465.rds")
ko_comm_a <- readRDS("ko_comm_a.rds")
ko_comm_b <- readRDS("ko_comm_b.rds") %>%
  dplyr::select(colnames(ko_comm_a))
ko_comm_t <- readRDS("ko_comm_t.rds")
ko_comm_o <- readRDS("ko_comm_o.rds")

sum(rownames(ko_comm_a) != f_meta_MUSiCC$taxon_oid)
sum(rownames(ko_comm_b) != f_meta_MUSiCC$taxon_oid)
sum(rownames(ko_comm_o) != f_meta_MUSiCC$taxon_oid)
sum(rownames(ko_comm_tot) != f_meta_MUSiCC$taxon_oid)

# Perform CPM assembled reads normalization
ko_comm_a_cpm <- ko_comm_a %>%
  mutate("K00194" = (K00194 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00399" = (K00399 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00440" = (K00440 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14080" = (K14080 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14084" = (K14084 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16177" = (K16177 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16179" = (K16179 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16954" = (K16954 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`)

ko_comm_b_cpm <- ko_comm_b %>%
  mutate("K00194" = (K00194 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00399" = (K00399 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00440" = (K00440 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14080" = (K14080 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14084" = (K14084 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16177" = (K16177 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16179" = (K16179 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16954" = (K16954 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`)

ko_comm_o_cpm <- ko_comm_o %>%
  mutate("K00194" = (K00194 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00399" = (K00399 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00440" = (K00440 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14080" = (K14080 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14084" = (K14084 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16177" = (K16177 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16179" = (K16179 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16954" = (K16954 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`)

ko_comm_t_cpm <- ko_comm_t %>%
  mutate("K00194" = (K00194 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00399" = (K00399 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K00440" = (K00440 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14080" = (K14080 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K14084" = (K14084 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16177" = (K16177 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16179" = (K16179 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`,
         "K16954" = (K16954 * 1000000)/f_meta_MUSiCC$`Genome Size   * assembled`)

# Data frames
gene_plot_a <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                          "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                          "Taxonomy" = "Archaea",
                          "mcrA" = ko_comm_a_cpm$K00399,
                          "cdhD" = ko_comm_a_cpm$K00194,
                          "frhA" = ko_comm_a_cpm$K00440,
                          "mttC" = ko_comm_a_cpm$K14084,
                          "mtbC" = ko_comm_a_cpm$K16179,
                          "mtmC" = ko_comm_a_cpm$K16177,
                          "mtaA" = ko_comm_a_cpm$K14080,
                          "mtsA" = ko_comm_a_cpm$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))

gene_plot_b <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                          "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                          "Taxonomy" = "Bacteria",
                          "mcrA" = ko_comm_b_cpm$K00399,
                          "cdhD" = ko_comm_b_cpm$K00194,
                          "frhA" = ko_comm_b_cpm$K00440,
                          "mttC" = ko_comm_b_cpm$K14084,
                          "mtbC" = ko_comm_b_cpm$K16179,
                          "mtmC" = ko_comm_b_cpm$K16177,
                          "mtaA" = ko_comm_b_cpm$K14080,
                          "mtsA" = ko_comm_b_cpm$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))

gene_plot_o <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                          "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                          "Taxonomy" = "Other",
                          "mcrA" = ko_comm_o_cpm$K00399,
                          "cdhD" = ko_comm_o_cpm$K00194,
                          "frhA" = ko_comm_o_cpm$K00440,
                          "mttC" = ko_comm_o_cpm$K14084,
                          "mtbC" = ko_comm_o_cpm$K16179,
                          "mtmC" = ko_comm_o_cpm$K16177,
                          "mtaA" = ko_comm_o_cpm$K14080,
                          "mtsA" = ko_comm_o_cpm$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))

gene_plot_t <- data.frame("Ecosystem" = f_meta_MUSiCC$Ecosystem,
                          "Hypothesis" = f_meta_MUSiCC$Hypothesis,
                          "Taxonomy" = "Total",
                          "mcrA" = ko_comm_t_cpm$K00399,
                          "cdhD" = ko_comm_t_cpm$K00194,
                          "frhA" = ko_comm_t_cpm$K00440,
                          "mttC" = ko_comm_t_cpm$K14084,
                          "mtbC" = ko_comm_t_cpm$K16179,
                          "mtmC" = ko_comm_t_cpm$K16177,
                          "mtaA" = ko_comm_t_cpm$K14080,
                          "mtsA" = ko_comm_t_cpm$K16954) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline"))

gene_plot_cpm <- rbind(gene_plot_a, gene_plot_b, gene_plot_o)



#### __ (III) Stats ####
# Run a loop for each
kruskal_results_genes_a <- as.data.frame(matrix(data = NA, 11, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 4:11) {
  k <- kruskal.test(gene_plot_a[[i]] ~ gene_plot_a$Ecosystem)
  kruskal_results_genes_a$Gene[i] <- names(gene_plot_a)[i]
  kruskal_results_genes_a$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes_a$P[i] <- k$p.value
}
kruskal_results_genes_a <- kruskal_results_genes_a %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"),
         Taxonomy = "Archaea") %>%
  select(Taxonomy, everything(), -P)

kruskal_results_genes_b <- as.data.frame(matrix(data = NA, 11, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 4:11) {
  k <- kruskal.test(gene_plot_b[[i]] ~ gene_plot_b$Ecosystem)
  kruskal_results_genes_b$Gene[i] <- names(gene_plot_b)[i]
  kruskal_results_genes_b$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes_b$P[i] <- k$p.value
}
kruskal_results_genes_b <- kruskal_results_genes_b %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"),
         Taxonomy = "Bacteria") %>%
  select(Taxonomy, everything(), -P)

kruskal_results_genes_o <- as.data.frame(matrix(data = NA, 11, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 4:11) {
  k <- kruskal.test(gene_plot_o[[i]] ~ gene_plot_o$Ecosystem)
  kruskal_results_genes_o$Gene[i] <- names(gene_plot_o)[i]
  kruskal_results_genes_o$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes_o$P[i] <- k$p.value
}
kruskal_results_genes_o <- kruskal_results_genes_o %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"),
         Taxonomy = "Other") %>%
  select(Taxonomy, everything(), -P)

kruskal_results_genes_t <- as.data.frame(matrix(data = NA, 11, 3)) %>%
  set_names(c("Gene", "X2", "P"))
for (i in 4:11) {
  k <- kruskal.test(gene_plot_t[[i]] ~ gene_plot_t$Ecosystem)
  kruskal_results_genes_t$Gene[i] <- names(gene_plot_t)[i]
  kruskal_results_genes_t$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_genes_t$P[i] <- k$p.value
}
kruskal_results_genes_t <- kruskal_results_genes_t %>%
  filter(is.na(Gene) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"),
         Taxonomy = "All") %>%
  select(Taxonomy, everything(), -P)

kruskal_results <- rbind(kruskal_results_genes_a,
                         kruskal_results_genes_b,
                         kruskal_results_genes_o,
                         kruskal_results_genes_t)

# Each ecosystem
long <- gene_plot_a %>%
  pivot_longer(., cols = c("mcrA","cdhD","frhA","mttC",
                           "mtbC","mtmC","mtaA","mtsA"))

landfill <- subset(long, Ecosystem == "Landfill")
leveneTest(value ~ name, data = landfill)
kruskal.test(value ~ name, data = landfill)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = landfill)

sew <- subset(long, Ecosystem == "Sewage treatment")
leveneTest(value ~ name, data = sew)
kruskal.test(value ~ name, data = sew)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = sew)

rice <- subset(long, Ecosystem == "Rice field")
leveneTest(value ~ name, data = rice) # Good
kruskal.test(value ~ name, data = rice)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = rice)

wetF <- subset(long, Ecosystem == "Wetland (freshwater)")
leveneTest(value ~ name, data = wetF)
kruskal.test(value ~ name, data = wetF)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = wetF)

oce <- subset(long, Ecosystem == "Ocean")
leveneTest(value ~ name, data = oce)
kruskal.test(value ~ name, data = oce)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = oce)

huli <- subset(long, Ecosystem == "Humans and livestock")
leveneTest(value ~ name, data = huli)
kruskal.test(value ~ name, data = huli)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = huli)

term <- subset(long, Ecosystem == "Termites")
leveneTest(value ~ name, data = term)
kruskal.test(value ~ name, data = term)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = term)

vent <- subset(long, Ecosystem == "Hydrothermal vent")
leveneTest(value ~ name, data = vent)
kruskal.test(value ~ name, data = vent)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = vent)

wetC <- subset(long, Ecosystem == "Wetland (coastal)")
leveneTest(value ~ name, data = wetC)
kruskal.test(value ~ name, data = wetC)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = wetC)

hsal <- subset(long, Ecosystem == "Hypersaline")
leveneTest(value ~ name, data = hsal)
kruskal.test(value ~ name, data = hsal)
kwAllPairsNemenyiTest(value ~ as.factor(name), data = hsal)



#### __ (IV) Graph ####
gene_plot_cpm_long <- gene_plot_cpm %>%
  pivot_longer(c("mcrA", "cdhD", "frhA", "mttC", 
                 "mtbC", "mtmC", "mtaA", "mtsA"), 
               names_to = "Gene", values_to = "Abundance") %>%
  mutate(Gene = as.factor(Gene)) %>%
  mutate(Gene = recode_factor(Gene,
                              "mcrA" = "mcrA",
                              "cdhD" = "cdhD",
                              "frhA" = "frhA",
                              "mttC" = "mttC",
                              "mtbC" = "mtbC", 
                              "mtmC" = "mtmC", 
                              "mtaA" = "mtaA",
                              "mtsA" = "mtsA")) %>%
  mutate(Pathway = recode_factor(Gene,
                                 "mcrA" = "All (mcrA)",
                                 "cdhD" = "Acetate (cdhD)",
                                 "frhA" = "H2/CO2 (frhA)",
                                 "mtaA" = "Methanol (mtaA)",
                                 "mttC" = "TMA (mttC)",
                                 "mtbC" = "DMA (mtbC)",
                                 "mtmC" = "MMA (mtmC)",
                                 "mtsA" = "DMS/MeSH/MMPA (mtsA)")) %>%
  droplevels() %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetland (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methylotrophic"))

lyu_only_summary_m <- gene_plot_cpm_long %>%
  group_by(Hypothesis, Ecosystem, Taxonomy, Gene, Pathway) %>%
  summarise(mean = mean(Abundance),
            se = std.error(Abundance))

# Facet
facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethyl-based",
                 "Archaea" = "Archaea", 
                 "Bacteria" = "Bacteria",
                 "Other" = "Other")
pdf("Figure4_CPMforppt_facet.pdf", width = 9, height = 5)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = position_dodge(0.75)) +
  geom_linerange(aes(x = Ecosystem, ymin = mean - se, ymax = mean + se, 
                     group = Pathway),
                 position = position_dodge(0.75)) +
  labs(x = NULL, 
       y = "CPM") +
  scale_fill_manual(values = c(viridis(20)[20],
                               viridis(20)[15],
                               viridis(20)[10],
                               viridis(20)[5],
                               brewer.pal(4, "Purples")[1],
                               brewer.pal(4, "Purples")[2],
                               brewer.pal(4, "Purples")[3],
                               brewer.pal(4, "Purples")[4])) +
  facet_grid(Taxonomy ~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))
dev.off()

# Nested facet
lyu_only_summary_m$Taxonomy <- factor(lyu_only_summary_m$Taxonomy,
                                      levels = c("Other", "Bacteria", "Archaea"))
lyu_only_summary_m <- lyu_only_summary_m %>%
  arrange(Taxonomy)
levels(lyu_only_summary_m$Ecosystem)
facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethyl-based",
                 "Landfill" = "Landfill", 
                 "Sewage treatment" = "Sewage\ntreatment",
                 "Rice field" = "Rice\nfield",
                 "Wetland (freshwater)" = "Wetland\n(freshwater)",
                 "Ocean" = "Ocean",
                 "Humans and livestock" = "Humans\nand\nlivestock",
                 "Termites" = "Termites",
                 "Hydrothermal vent" = "Hydrothermal\nvent",
                 "Wetland (coastal)" = "Wetland\n(coastal)",
                 "Hypersaline" = "Hypersaline")
pdf("Figure4_CPMforppt_facetnested.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Pathway, mean, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "CPM") +
  facet_nested(~ Hypothesis + Ecosystem, scales = "free_x", space = "free",
               labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.00, 0.01)) +  
  scale_fill_manual(values = c("grey", "blue", "red")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                   margin = margin(0,0,0,0, "pt")),
        strip.text = element_text(size = 5.5),
        strip.background = element_rect(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.1, 0.7),
        legend.key.size = unit(0.3, "cm"))
dev.off()

pdf("Figure4_CPMforppt_facetnested_noleg.pdf", width = 9, height = 4)
ggplot(lyu_only_summary_m, aes(Pathway, mean, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "CPM") +
  facet_nested(~ Hypothesis + Ecosystem, scales = "free_x", space = "free",
               labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.00, 0.00)) +  
  scale_fill_manual(values = c("grey", "blue", "red")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                   margin = margin(0,0,0,0, "pt")),
        strip.text = element_text(size = 5.5),
        strip.background = element_rect(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.position = "none")
dev.off()

# Each gene
cdhd <- subset(lyu_only_summary_m, Gene == "cdhD")
ggplot(cdhd, aes(Ecosystem, mean, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "CPM") +
  facet_nested(~ Hypothesis + Ecosystem, scales = "free_x", space = "free",
               labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.00, 0.01)) +  
  scale_fill_manual(values = c("grey", "blue", "red")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                   margin = margin(0,0,0,0, "pt")),
        strip.text = element_text(size = 5.5),
        strip.background = element_rect(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.1, 0.7),
        legend.key.size = unit(0.3, "cm"))

frhA <- subset(lyu_only_summary_m, Gene == "frhA")
ggplot(frhA, aes(Ecosystem, mean, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "CPM") +
  facet_nested(~ Hypothesis + Ecosystem, scales = "free_x", space = "free",
               labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.00, 0.01)) +  
  scale_fill_manual(values = c("grey", "blue", "red")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                   margin = margin(0,0,0,0, "pt")),
        strip.text = element_text(size = 5.5),
        strip.background = element_rect(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.1, 0.7),
        legend.key.size = unit(0.3, "cm"))


mb <- subset(lyu_only_summary_m, 
             Gene == "mttC" |
               Gene == "mtbC" |
               Gene == "mtmC" |
               Gene == "mtaA" |
               Gene == "mtsA")
ggplot(mb, aes(Pathway, mean, fill = Taxonomy)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, 
       y = "CPM") +
  facet_nested(~ Hypothesis + Ecosystem, scales = "free_x", space = "free",
               labeller = as_labeller(facet_names)) +
  scale_y_continuous(expand = c(0.00, 0.01)) +  
  scale_fill_manual(values = c("grey", "blue", "red")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                   margin = margin(0,0,0,0, "pt")),
        strip.text = element_text(size = 5.5),
        strip.background = element_rect(size = 0.2),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.position = c(0.1, 0.7),
        legend.key.size = unit(0.3, "cm"))


# Stacked (difficult, not really working/accurate!)
facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methylotrophic" = "Hypothesized\nMethyl-based")
lyu_only_summary_m$Taxonomy <- factor(lyu_only_summary_m$Taxonomy,
                                      levels = c("Other", "Bacteria", "Archaea"))
lyu_only_summary_m <- lyu_only_summary_m %>%
  arrange(Taxonomy)
ggplot(lyu_only_summary_m, aes(Ecosystem, mean, fill = Pathway, group = Pathway)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.1) +
  labs(x = NULL, 
       y = "CPM") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))

ggplot(lyu_only_summary_m, 
       aes(x = interaction(Pathway,Ecosystem), y = mean, fill = Taxonomy)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  labs(x = NULL, 
       y = "CPM") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.box.margin=margin(0,-5,0,60))



#### (b) Taxonomic (family) ####
# Metadata for mctools import
metad <- f_meta_MUSiCC %>%
  select(sampleID, everything())
# write.table(metad, file = "~/Desktop/Review/f_meta_MUSiCC.txt", sep = "\t", row.names = F)

# Import
tax_table_fp <- file.path("~/Desktop/Review/IMG_families_mctoolsr.txt")
map_fp <- file.path("~/Desktop/Review/f_meta_MUSiCC.txt")
input = load_taxa_table(tax_table_fp, map_fp)
input$map_loaded <- input$map_loaded %>%
  mutate(sampleID = paste("X", taxon_oid, sep = ""))

# Match the functional analysis (should already be fine)
input = filter_data(input,
                    filter_cat = "sampleID",
                    keep_vals = f_meta_MUSiCC$sampleID)

# Check sequencing depth 
sort(colSums(input$data_loaded))
mean(colSums(input$data_loaded))
se(colSums(input$data_loaded))
input$map_loaded$count <- colSums(input$data_loaded)
ggplot(input$map_loaded, aes(reorder(Ecosystem, count, mean), count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "# Reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

ggplot(input$map_loaded, aes(`Genome Size   * assembled`, count)) +
  geom_point(size = 1.5, alpha = 0.25) +
  labs(x = "Assembled genome size", 
       y = "Assigned family reads") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10))

# Rename size column
input$map_loaded$GenomeSize = input$map_loaded$`Genome Size   * assembled`

# Phyla
tax_sum_phyla <- summarize_taxonomy(input, level = 2, report_higher_tax = T, relative = FALSE)
# Note there are no Bathyarchaeota or Verstraetearchaeota, but check class too
# Class
tax_sum_class <- summarize_taxonomy(input, level = 3, report_higher_tax = T, relative = FALSE)
# No Bathyarchaeia

# Summarize by family, extract methanogens, calculate CPM
tax_sum_family_wTax <- summarize_taxonomy(input, level = 5, report_higher_tax = T, relative = FALSE)
methano_wTax <- tax_sum_family_wTax[grep("Methano", rownames(tax_sum_family_wTax)),]
methano_wTax <- methano_wTax[!grepl("Plasmid", rownames(methano_wTax)),]
# Note this also contains a methanotroph Methermicoccaceae (also methanogen?), not captured below though
# Careful with Methanoperedens too

# Without higher tax and no plasmid
tax_sum_family <- summarize_taxonomy(input, level = 5, report_higher_tax = FALSE, relative = FALSE)
methano <- tax_sum_family[grep("Methano", rownames(tax_sum_family)),]
tax_sum_family_df <- t(methano) %>%
  as.data.frame() %>%
  mutate(`Candidatus Methanoperedenaceae` = (`Candidatus Methanoperedenaceae`*1000000)/input$map_loaded$GenomeSize,
         Methanobacteriaceae = (Methanobacteriaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocaldococcaceae = (Methanocaldococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocellaceae = (Methanocellaceae*1000000)/input$map_loaded$GenomeSize,
         Methanococcaceae = (Methanococcaceae*1000000)/input$map_loaded$GenomeSize,
         Methanocorpusculaceae = (Methanocorpusculaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomassiliicoccaceae = (Methanomassiliicoccaceae*1000000)/input$map_loaded$GenomeSize,
         Methanomicrobiaceae = (Methanomicrobiaceae*1000000)/input$map_loaded$GenomeSize,
         Methanonatronarchaeaceae = (Methanonatronarchaeaceae*1000000)/input$map_loaded$GenomeSize,
         Methanopyraceae = (Methanopyraceae*1000000)/input$map_loaded$GenomeSize,
         Methanoregulaceae = (Methanoregulaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosaetaceae = (Methanosaetaceae*1000000)/input$map_loaded$GenomeSize,
         Methanosarcinaceae = (Methanosarcinaceae*1000000)/input$map_loaded$GenomeSize,
         Methanospirillaceae = (Methanospirillaceae*1000000)/input$map_loaded$GenomeSize,
         Methanothermaceae = (Methanothermaceae*1000000)/input$map_loaded$GenomeSize,
         Methanotrichaceae = (Methanotrichaceae*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanomicrobiales` = (`unclassified Methanomicrobiales`*1000000)/input$map_loaded$GenomeSize,
         `unclassified Methanosarcinales` = (`unclassified Methanosarcinales`*1000000)/input$map_loaded$GenomeSize) 

input$map_loaded <- input$map_loaded %>%
  mutate(Ecosystem = as.factor(Ecosystem),
         Hypothesis = as.factor(Hypothesis)) %>%
  mutate(Ecosystem = recode_factor(Ecosystem,
                                   "Landfill" = "Landfill",
                                   "Sewage treatment" = "Sewage treatment",
                                   "Rice field" = "Rice field",
                                   "Wetland (freshwater)" = "Wetlands (freshwater)",
                                   "Ocean" = "Ocean",
                                   "Cow gut" = "Humans and livestock",
                                   "Human gut" = "Humans and livestock",
                                   "Sheep gut" = "Humans and livestock",
                                   "Termite gut" = "Termites",
                                   "Hydrothermal vent" = "Hydrothermal vent",
                                   "Wetland (coastal)" = "Wetland (coastal)",
                                   "Hypersaline" = "Hypersaline")) %>%
  mutate(Hypothesis = recode_factor(Hypothesis,
                                    "H2, Ac" = "Acetoclastic",
                                    "H2" = "Hydrogenotrophic",
                                    "Me" = "Methyl-based"))
sum(input$map$sampleID != rownames(tax_sum_family_df))
input$map_loaded <- cbind(input$map_loaded, tax_sum_family_df)



#### __(I) Stats ####
# Run a loop
kruskal_results_taxa <- as.data.frame(matrix(data = NA, 44, 3)) %>%
  set_names(c("Family", "X2", "P"))
for (i in 27:44) {
  k <- kruskal.test(input$map_loaded[[i]] ~ input$map_loaded$Ecosystem)
  kruskal_results_taxa$Family[i] <- names(input$map_loaded)[i]
  kruskal_results_taxa$X2[i] <- round(k$statistic, digits = 2)
  kruskal_results_taxa$P[i] <- k$p.value
}
kruskal_results_taxa <- kruskal_results_taxa %>%
  filter(is.na(Family) == F) %>%
  mutate(Pfdr = p.adjust(P, method = "fdr"))


#### __(II) Graph ####
family_plot_long <- pivot_longer(input$map_loaded,
                                 cols = names(input$map_loaded)[27:44],
                                 names_to = "Family",
                                 values_to = "Abundance")
test <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  group_by(Family) %>%
  summarise(max_abund = max(mean_abund))

test2 <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  group_by(Family) %>%
  summarise(min_abund = min(mean_abund))

test3 <- family_plot_long %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(Abundance))

test4 <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance))

bars <- family_plot_long %>%
  group_by(Ecosystem, Hypothesis, Family) %>%
  summarise(mean_abund = mean(Abundance)) %>%
  mutate(Family = as.factor(Family)) %>%
  filter(Family != "unclassified Methanomicrobiales",
         Family != "unclassified Methanosarcinales",
         Family != "Candidatus Methanoperedenaceae")

filter_bars <- bars %>%
  group_by(Family) %>%
  summarise(max_abund = max(mean_abund)) %>%
  filter(max_abund > 1)

bars_top <- bars %>%
  filter(Family %in% filter_bars$Family)

sort_bars <- bars_top %>%
  group_by(Family) %>%
  summarise(sum_abund = sum(mean_abund)) %>%
  arrange(sum_abund)

bars_top_sorted <- bars_top %>%
  mutate(Family = factor(Family,
                         levels = sort_bars$Family))

bars_top_sorted_guild <- bars_top_sorted %>%
  mutate(Family = factor(Family,
                         levels = c("Methanosaetaceae", "Methanotrichaceae",
                                    "Methanospirillaceae", "Methanocellaceae",
                                    "Methanocorpusculaceae", "Methanococcaceae",
                                    "Methanomicrobiaceae", "Methanocaldococcaceae",
                                    "Methanobacteriaceae", "Methanoregulaceae",
                                    "Methanomassiliicoccaceae", "Methanosarcinaceae")))

n_families <- bars %>%
  group_by(Hypothesis, Ecosystem) %>%
  summarise(n_fam = sum(mean_abund > 0),
            tot = sum(mean_abund)) %>%
  mutate(y = tot + 2)

nb.cols <- nrow(filter_bars)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

facet_names <- c("Acetoclastic" = "Hypothesized\nAcetoclastic",
                 "Hydrogenotrophic" = "Hypothesized\nHydrogenotrophic",
                 "Methyl-based" = "Hypothesized\nMethyl-based")

# Figure 4 - finalize in powerpoint
pdf("Figure5_forppt.pdf", width = 8.5, height = 4.5)
ggplot(bars_top_sorted_guild, aes(Ecosystem, mean_abund, fill = Family)) +
  geom_bar(stat = "identity", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.box.margin=margin(0,-5,0,15))
dev.off()

# Remake figure aggregated by guild
bars_top_sorted_guild$Guild <- NA
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanosaetaceae") {
    bars_top_sorted_guild$Guild[i] <- "Acetoclastic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanotrichaceae") {
    bars_top_sorted_guild$Guild[i] <- "Acetoclastic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanospirillaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocellaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocorpusculaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanococcaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanomicrobiaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanocaldococcaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanobacteriaceae") {
    bars_top_sorted_guild$Guild[i] <- "Mixotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanoregulaceae") {
    bars_top_sorted_guild$Guild[i] <- "Hydrogenotrophic"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanomassiliicoccaceae") {
    bars_top_sorted_guild$Guild[i] <- "Methyl-reducing"
  }
}
for (i in 1:nrow(bars_top_sorted_guild)) {
  if (bars_top_sorted_guild$Family[i] == "Methanosarcinaceae") {
    bars_top_sorted_guild$Guild[i] <- "Mixotrophic"
  }
}
bars_top_sorted_guild$Guild <- as.factor(bars_top_sorted_guild$Guild)
bars_top_sorted_guild_ag <- bars_top_sorted_guild %>%
  group_by(Ecosystem, Hypothesis, Guild) %>%
  summarise(mean_abund = sum(mean_abund))


# Guild figure
pdf("Manuscript/FigureS3.pdf", width = 8.5, height = 4.5)
ggplot(bars_top_sorted_guild_ag, aes(Ecosystem, mean_abund, fill = Guild)) +
  geom_bar(stat = "identity", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7))
dev.off()

png("Manuscript/FigureS3.png", width = 8.5, height = 4.5, unit = "in", res = 300)
ggplot(bars_top_sorted_guild_ag, aes(Ecosystem, mean_abund, fill = Guild)) +
  geom_bar(stat = "identity", size = 0.25) +
  geom_text(data = n_families,
            aes(x = Ecosystem, y = y, label = n_fam),
            inherit.aes = F, size = 3) +
  labs(x = NULL, 
       y = "Mean abundance (CPM)") +
  facet_grid(~ Hypothesis, scales = "free_x", space = "free", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 7))
dev.off()

hue_pal()(4)

# Sample size
table(input$map_loaded$Ecosystem)

# Check genera
tax_sum_genera <- summarize_taxonomy(input, level = 6, report_higher_tax = T, relative = FALSE)
# Oh yea, this was a family level download, so no genera...


#### Map ####
# Sample map with ggplot
world <- map_data("world")

pdf("Manuscript/FigureS1.pdf", width = 7, height = 5)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  geom_point(data = f_meta_MUSiCC, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

png("Manuscript/FigureS1.png", width = 7, height = 5, units = "in", res = 300)
ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "lightgray", size = 0.1) +
  geom_point(data = f_meta, aes(x = Longitude, y = Latitude),
             size = 1, color = "red", alpha = 0.5) +
  theme_void() +
  labs(x = NULL,
       y = NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



#### Table S1 ####
# Save Updated Table S1 with the 395 samples actually analyzed
# -Combine and remove unnecessary columns
# -Add additional data status column from GOLD
# -Reorder columns
# Add Permission notes from the previous version (with n = 395)
tableS1_notes <- tableS1 %>%
  select(taxon_oid, PermissionNotes)

metadata <- f_meta_MUSiCC %>%
  select(-sampleID) %>%
  select(-`Contact Name`, -`Contact Email`, -`JGI Data Utilization Status`, -`JGI Data Utilization Status_GOLD_07Apr2022`, -PermissionNotes, everything(), `Contact Name`, `Contact Email`, `JGI Data Utilization Status`, `JGI Data Utilization Status_GOLD_07Apr2022`) %>%
  select(-PermissionNotes) %>%
  left_join(., tableS1_notes, by = "taxon_oid")

# Export as Supplementary Table 1
# write_xlsx(metadata, "TableS1_Updated_465.xlsx", format_headers = F)

# In Excel fill PermissionNotes column confirming instances in which a PI was contacted and gave permission

# Get the extra 70 that were added
t395 <- read_xlsx("~/Desktop/Review/TableS1_Updated_395.xlsx")
t465 <- read_xlsx("~/Desktop/Review/Manuscript/TableS1_Updated_465.xlsx")
added_oid <- as.data.frame(t465$taxon_oid) %>%
  set_names("taxonoid") %>%
  filter(taxonoid %notin% t395$taxon_oid)
write.csv(added_oid, "added70.csv")