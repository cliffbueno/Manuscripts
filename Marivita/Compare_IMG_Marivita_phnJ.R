# Marivita and phnJ comparative analysis
# Streamlined version of Compare_IMG_736.R analysis just for this paper
# In the end produces Figure 4 for manuscript

#### Retrieving Data ####
# To retrieve data from IMG, first make a genome set with all of the IMG ID's from the literature or that are from saline environments. Add features to the table to download metadata.
# To retrieve functional data (KO), use the statistical analysis tools.
# To retrieve taxonomic data (class, family, genus), use the statistical analysis tools.
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
suppressWarnings(suppressMessages(library(PMCMR))) # For Nemenyi posthoc test
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

setwd("~/Desktop/MetagenomeComp")
options(max.print = 2000000)
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullj <- function(df) df[chull(df$Axis01j, df$Axis02j),]

# Functional table (KO)
f <- read.delim("stat_ko736_18-jan-2021/UI_data_output.txt") %>%
  separate(Feature, into = c("Text", "KO"), sep = ":") %>%
  column_to_rownames(var = "KO") %>%
  select(745:1480)
names(f) <- abbreviate(names(f), minlength = 11)
f[is.na(f)] <- 0
f_comm <- as.data.frame(t(f)) %>%
  arrange(row.names(.))

# Metadata
f_meta <- read.delim("SalineSoilSedMetagenomes_IMG_736.txt") %>%
  arrange(sampleID) %>%
  mutate_if(is.character, as.factor)
sum(f_meta$sampleID != rownames(f_comm))
f_meta$richness_KO = specnumber(f_comm)
min(f_meta$richness_KO) # 1 zero, remove.
f_comm <- f_comm %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(richness_KO = f_meta$richness_KO) %>%
  filter(richness_KO > 0) %>%
  select(-richness_KO)
f_meta <- subset(f_meta, richness_KO > 0)
sum(f_meta$sampleID != f_comm$sampleID)
rownames(f_comm) <- f_comm$sampleID
f_comm$sampleID <- NULL

# Remove water and rice paddy samples. Also Mud flat and mud volcano (only 1 sample each)
f_comm <- f_comm %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(EcosystemL3 = f_meta$EcosystemL3) %>%
  filter(EcosystemL3 != "Delete") %>%
  filter(EcosystemL3 != "Mud_flat") %>%
  filter(EcosystemL3 != "Mud_volcano") %>%
  select(-EcosystemL3)
f_meta <- f_meta %>%
  filter(EcosystemL3 != "Delete") %>%
  filter(EcosystemL3 != "Mud_flat") %>%
  filter(EcosystemL3 != "Mud_volcano")
sum(f_meta$sampleID != f_comm$sampleID)
rownames(f_comm) <- f_comm$sampleID
f_comm$sampleID <- NULL
f_comm$EcosystemL3 <- NULL
f_meta <- droplevels(f_meta)
saveRDS(f_meta, "~/Desktop/Marivita/f_meta.rds")

# DESeq
dds <- DESeqDataSetFromMatrix(countData = t(f_comm) + 1,
                              colData = f_meta,
                              design = ~ Study.Name)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
f_comm_DESeq <- as.data.frame(t(counts(dds, normalized = T)))
saveRDS(f_comm_DESeq, "~/Desktop/Marivita/phnJ_DESeq.rds")

#### Start Here ####
f_comm_DESeq <-readRDS("~/Desktop/Marivita/phnJ_DESeq.rds")
f_meta <- readRDS("~/Desktop/Marivita/f_meta.rds")


#### phnJ ####

# Check IDs
sum(f_meta$sampleID != rownames(f_comm_DESeq))

# Plot mcrA and phnJ by ecosystem
gene_plot <- data.frame("EcosystemL3" = f_meta$EcosystemL3,
                        "log2DESeqmcrA" = log2(f_comm_DESeq$K00399),
                        "log2DESeqphnJ" = log2(f_comm_DESeq$K06163),
                        "Database" = "IMG")

f_meta <- f_meta %>%
  mutate("log2DESeqmcrA" = log2(f_comm_DESeq$K00399),
         "log2DESeqphnJ" = log2(f_comm_DESeq$K06163))


leveneTest(log2DESeqphnJ ~ EcosystemL3, data = f_meta)
m <- aov(log2DESeqphnJ ~ EcosystemL3, data = f_meta)
summary(m)
shapiro.test(m$residuals)
TukeyHSD(m)
kruskal.test(log2DESeqphnJ ~ EcosystemL3, data = f_meta)
kwAllPairsNemenyiTest(log2DESeqphnJ ~ EcosystemL3, data = f_meta)

ggplot(f_meta, aes(reorder(EcosystemL4, log2DESeqphnJ, mean), log2DESeqphnJ)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "phnJ abundance (log2 DESeq count)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))

ggplot(f_meta, aes(reorder(EcosystemL4, log2DESeqmcrA, mean), log2DESeqmcrA)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "mcrA abundance (log2 DESeq count)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        plot.margin = margin(2, 2, 2, 15))

# Multipanel mcrA, phnJ for presentations 
gene_long <- melt(f_meta, 
                  id.vars = c("EcosystemL4"), 
                  measure.vars = c("log2DESeqmcrA", "log2DESeqphnJ"))
  
gene_long$EcosystemL4 <- factor(gene_long$EcosystemL4,
                                levels = c("Hot_spring", "Estuary", "Methane_seep",
                                        "Ocean", "Salt_marsh", "Restored_saltern", "Mangrove",
                                        "Ocean_coastal", "Hydrothermal_vent", "Coastal_wetland_SF",
                                        "Lagoon", "Coastal_wetland_R2A", "Coastal_wetland",
                                        "Hypersaline_lake", "Unrestored_saltern_R2",
                                        "Unrestored_saltern_R1"))
gene_long$variable <- factor(gene_long$variable,
                             levels = c("log2DESeqmcrA", "log2DESeqphnJ"))

facet_names <-  c(`log2DESeqmcrA` = "(a) mcrA (log2 DESeq2 normalized counts)",
                  `log2DESeqphnJ` = "(b) phnJ (log2 DESeq2 normalized counts)")

my_labeller <- as_labeller(c(log2DESeqmcrA = "(a)~mcrA~(log[2]~DESeq2~normalized~counts)", 
                             log2DESeqphnJ = "(b)~phnJ~(log[2]~DESeq2~normalized~counts)"),
                           default = label_parsed)

x_label <- c("Hot spring", "Estuary", "Methane seep", "Open ocean", "Salt marsh", "Restored saltern",
             "Mangrove", "Coastal ocean", "Hydrothermal vent", "Coastal wetland (SF)", "Lagoon",
             "Reference wetland", "Coastal wetland", "Hypersaline lake", "Unrestored saltern R2",
             "Unrestored saltern R1")               

pdf("~/Desktop/SouthBay/IMG_Compare_mcrA_phnJ.pdf", width = 7, height = 5)
ggplot(gene_long, aes(EcosystemL4, value, color = EcosystemL4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "Abundance") +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "red",
                                "black", "black", "black", "black", "black", "red",
                                "black", "black", "red", "red")) +
  scale_x_discrete(labels = x_label) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", labeller = my_labeller) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(size = 10))
dev.off()


#### Marivita ####
# Import
tax_table_fp <- file.path("~/Desktop/MetagenomeComp/IMG_otutable_wTax_mctoolsr_736.txt")
map_fp <- file.path("~/Desktop/MetagenomeComp/SalineSoilSedMetagenomes_IMG_736.txt")
input = load_taxa_table(tax_table_fp, map_fp)

sum(is.na(input$map_loaded$`Genome Size   * assembled`))

# Filter out non-sediment (water) samples
input = filter_data(input, "EcosystemL3", "Delete")

# For all data get rid of samples with less than 1000 genus reads
count <- as.data.frame(sort(colSums(input$data_loaded))) %>%
  filter(`sort(colSums(input$data_loaded))` < 1000)
input <- filter_data(input,
                     filter_cat = "sampleID",
                     filter_vals = rownames(count))

# Get rid of water samples
input <- filter_data(input,
                     filter_cat = "ClusterL1",
                     filter_vals = "Water")

# Get rid of ecosystems with 1 sample
count_ecosystem <- as.data.frame(table(input$map_loaded$EcosystemL1)) %>%
  filter(Freq == 1)
input <- filter_data(input,
                     filter_cat = "EcosystemL1",
                     filter_vals = count_ecosystem$Var1) # 590 samples

# Extract Marivita counts and transform to CPM
tax_sum_genus <- summarize_taxonomy(input, level = 6, report_higher_tax = FALSE)
tax_sum_genus <- summarize_taxonomy(input, level = 6, report_higher_tax = FALSE, relative = FALSE)
tax_sum_genus_df <- t(tax_sum_genus) %>%
  as.data.frame() %>%
  select(Marivita) %>%
  mutate(CPM = (Marivita * 1000000)/input$map_loaded$`Genome Size   * assembled`)
sum(rownames(tax_sum_genus_df) != rownames(input$map_loaded))

input$map_loaded$Marivita <- tax_sum_genus_df$CPM
leveneTest(input$map_loaded$Marivita ~ input$map_loaded$EcosystemL4)
m <- aov(input$map_loaded$Marivita ~ input$map_loaded$EcosystemL4)
summary(m)
shapiro.test(m$residuals)
TukeyHSD(m)
kruskal.test(input$map_loaded$Marivita ~ input$map_loaded$EcosystemL4)
kwAllPairsNemenyiTest(input$map_loaded$Marivita ~ input$map_loaded$EcosystemL4)
ggplot(input$map_loaded, aes(reorder(EcosystemL4, Marivita, mean), log2(Marivita+0.0000000001))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "Marivita (Log CPM Assembled reads)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"))



#### Multipanel Figure ####
f_meta$Variable <- "phnJ"
f_meta$Value <- f_meta$log2DESeqphnJ
y <- select(f_meta, Variable, EcosystemL4, Value)
input$map_loaded$Value <- log2(input$map_loaded$Marivita+0.0000000001)
input$map_loaded$Variable <- "Marivita"
input$map_loaded$EcosystemL4 <- as.factor(input$map_loaded$EcosystemL4)
z <- select(input$map_loaded, Variable, EcosystemL4, Value)

# Stats
leveneTest(input$map_loaded$Value ~ input$map_loaded$EcosystemL4)
m <- aov(input$map_loaded$Value ~ input$map_loaded$EcosystemL4)
summary(m)
shapiro.test(m$residuals)
TukeyHSD(m)
kruskal.test(input$map_loaded$Value ~ input$map_loaded$EcosystemL4)
nem <- kwAllPairsNemenyiTest(input$map_loaded$Value ~ input$map_loaded$EcosystemL4)

# Graph
d_long <- rbind(y, z)
d_long$EcosystemL4 <- factor(d_long$EcosystemL4,
                             levels = c("Hot_spring", "Estuary", "Methane_seep",
                             "Ocean", "Salt_marsh", "Restored_saltern", "Mangrove",
                             "Ocean_coastal", "Hydrothermal_vent", "Coastal_wetland_SF",
                             "Lagoon", "Coastal_wetland_R2A", "Coastal_wetland",
                             "Hypersaline_lake", "Unrestored_saltern_R2", "Unrestored_saltern_R1"))
d_long$Variable <- factor(d_long$Variable,
                          levels = c("phnJ", "Marivita"))

facet_names <-  c(`phnJ` = "(a) phnJ (log2 DESeq2 normalized counts)",
                  `Marivita` = "(b) Marivita (log2 counts per million)")

my_labeller <- as_labeller(c(phnJ = "(a)~phnJ~(log[2]~DESeq2~normalized~counts)", 
                             Marivita = "(b)~Marivita~(log[2]~counts~per~million)"),
                           default = label_parsed)

x_label <- c("Hot spring", "Estuary", "Methane seep", "Open ocean", "Salt marsh", "Restored saltern",
             "Mangrove", "Coastal ocean", "Hydrothermal vent", "Coastal wetland (SF)", "Lagoon",
             "Reference wetland", "Coastal wetland", "Hypersaline lake", "Unrestored saltern R2",
             "Unrestored saltern R1")

pdf("~/Desktop/Marivita/Manuscript/Figure4.pdf", width = 7, height = 5)
ggplot(d_long, aes(EcosystemL4, Value, color = EcosystemL4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "Abundance") +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "red",
                               "black", "black", "black", "black", "black", "red",
                               "black", "black", "red", "red")) +
  scale_x_discrete(labels = x_label) +
  facet_wrap(~Variable, ncol = 1, scales = "free_y", labeller = my_labeller) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(size = 10))
dev.off()

png("~/Desktop/Marivita/Manuscript/Figure4.png", width = 7, height = 5, unit = "in", res = 300)
ggplot(d_long, aes(EcosystemL4, Value, color = EcosystemL4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1.5, alpha = 0.25, width = 0.25) +
  labs(x = "Ecosystem", 
       y = "Abundance") +
  scale_color_manual(values = c("black", "black", "black", "black", "black", "red",
                                "black", "black", "black", "black", "black", "red",
                                "black", "black", "red", "red")) +
  scale_x_discrete(labels= x_label) +
  facet_wrap(~Variable, ncol = 1, scales = "free_y", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(size = 10))
dev.off()

# Quick check of the correlation (r = 0.36, p < 0.001)
c <- input$map_loaded %>%
  select(sampleID, Value) %>%
  left_join(., f_meta, by = "sampleID")
cor.test(c$Value.x, c$log2DESeqphnJ)
ggplot(c, aes(Value.x, log2DESeqphnJ)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Marivita",
       y = "phnJ") +
  theme_bw()

# Make Supplementary Table S1
tables1 <- f_meta %>%
  select(taxon_oid, EcosystemL4, Study.Name, Genome.Name...Sample.Name, Assembly.Method,
         Geographic.Location, Latitude, Longitude, Genome.Size.....assembled) %>%
  set_colnames(., c("IMG_ID", "Ecosystem", "Study_Name", "Genome_Name", "Assembly_Method",
                    "Geographic_Location", "Latitude", "Longitude", "Assembled_Genome_Size"))
# write_xlsx(tables1, path = "~/Desktop/Marivita/Manuscript/SupplementaryTable1.xlsx",
#           format_headers = FALSE)
