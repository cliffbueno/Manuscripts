# Heatmaps for KOs involved in maranogenesis and salt tolerance
# Genome level - 8 genomes
# - 2 MAGs and 6 other Marivita genomes
# Used BLAST Koala for KO assignment
# by Cliff Bueno de Mesquita, JGI, Spring 2021



#### Setup ####
library(plyr)
library(tidyverse)
library(readxl)
library(gplots)
library(scales)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(vegan)
library(RColorBrewer)
`%notin%` <- Negate(`%in%`)
show_col(hue_pal()(3))
show_col(viridis_pal()(3))
hue_pal()(13)[1]
save_pheatmap_pdf <- function(x, filename, width = 7, height = 5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

setwd("~/Desktop/Genomes/Marivita_KO_profiles_IMG_KOALA/")

# KO table - import 8 from BLAST Koala
MAG28_KO <- read.delim2("MarivitaKOs/Marivita28_82_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MAG28 = n()) %>%
  set_colnames(c("KO", "M. sp. SBSPR1"))

MAG10_KO <- read.delim2("MarivitaKOs/Marivita10_192_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MAG10 = n()) %>%
  set_colnames(c("KO", "M. sp. SBSPR2"))

MarCry_KO <- read.delim2("MarivitaKOs/MarCry_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MarCry = n()) %>%
  set_colnames(c("KO", "M. cryptomonadis"))

MarGeo_KO <- read.delim2("MarivitaKOs/MarGeo_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MarGeo = n()) %>%
  set_colnames(c("KO", "M. geojedonensis"))

MarHal_KO <- read.delim2("MarivitaKOs/MarHal_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MarHal = n()) %>%
  set_colnames(c("KO", "M. hallyeonensis"))

MarLac_KO <- read.delim2("MarivitaKOs/MarLac_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MarLac = n()) %>%
  set_colnames(c("KO", "M. lacus"))

MarLZ_KO <- read.delim2("MarivitaKOs/MarLZ_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MarLZ = n()) %>%
  set_colnames(c("KO", "M. sp. LZ-15-2"))

# Merge into master table of 1s and 0s
KO_table <- full_join(MAG10_KO, MAG28_KO, by = "KO") %>%
  full_join(., MarCry_KO, by = "KO") %>%
  full_join(., MarGeo_KO, by = "KO") %>%
  full_join(., MarHal_KO, by = "KO") %>%
  full_join(., MarLac_KO, by = "KO") %>%
  full_join(., MarLZ_KO, by = "KO") %>%
  column_to_rownames(var = "KO") %>%
  mutate_if(is.character, as.integer)
KO_table[is.na(KO_table)] <- 0
KO_table[KO_table > 0] <- 1
KO_table <- KO_table %>%
  rownames_to_column(var = "KO")
names(KO_table)[4] <- "M. cryptomonadis MP20-4"
names(KO_table)[8] <- "M. cryptomonadis LZ-15-2"



#### Methylphosphonate Degradation ####
# Note, all zeros in paths II and II so just show path I
a <- read_excel("~/Desktop/Marivita/MethylphosphonateKOs.xlsx", sheet = 1) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order)

mar_KOs <- left_join(a, KO_table, by = c("KEGG_KO" = "KO")) %>%
  replace(is.na(.), 0) %>%
  select(8:15)
names(mar_KOs)[4] <- "M. cryptomonadis MP20-4"
names(mar_KOs)[8] <- "M. cryptomonadis LZ-15-2"

mar_meta <- left_join(a, KO_table, by = c("KEGG_KO" = "KO")) %>%
  select(1:8)

mar_mat <- mar_KOs %>%
  column_to_rownames(var = "KO_def") %>%
  select(7, 3, 6, 2, 1, 4, 5) %>%
  as.matrix()

# Pretty heatmap
ann_rows <- data.frame(row.names = rownames(mar_mat), Pathway = mar_meta$Pathway_Specific)
ann_colors <- list(Pathway = c(`MPn synthesis` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[1],
                               `MPn transport` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[2],
                               `C-P Lyase` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[3],
                               `a-D-ribose 1,5-biphosphate` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[4],
                               `D-ribofuranose 5-phosphate` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[5],
                               `Non C-P Lyase` = hue_pal()(length(levels(mar_meta$Pathway_Specific)))[6]))
setwd("~/Desktop/Marivita/Manuscript/")
phm1 <- pheatmap(mar_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 fontsize = 10,
                 fontsize_row = 10,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 cluster_rows = F,
                 cluster_cols = F,
                 gaps_row = c(1, 4, 10, 11, 12))
save_pheatmap_pdf(phm1, "MethylphosphonateModulesKO_genome.pdf")



#### Salt Tolerance ####
salt_orig <- read_excel("~/Desktop/Methanolobus/SaltGenes.xlsx") %>%
  filter(KO != "NA") %>%
  mutate(Type = as.factor(Type),
         Solute = as.factor(Solute)) %>%
  group_by(KO) %>%
  slice_head(n = 1)
s_meta <- KO_table %>%
  filter(KO %in% salt_orig$KO) %>%
  left_join(., salt_orig, by = "KO") %>%
  arrange(Order) %>%
  droplevels()
salt <- read_excel("~/Desktop/Methanolobus/SaltGenes.xlsx") %>%
  filter(KO != "NA") %>%
  group_by(KO) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Order) %>%
  select(KO, Pathway_general, Code, Order)
s <- KO_table %>%
  filter(KO %in% salt$KO) %>%
  left_join(., salt, by = "KO") %>%
  mutate(KO = paste(KO, Code, sep = " ")) %>%
  mutate(KO = paste(KO, Pathway_general, sep = "; ")) %>%
  arrange(Order)
s_mat <- s %>%
  column_to_rownames(var = "KO") %>%
  select(7, 3, 6, 2, 1, 4, 5) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(s_mat), Type = s_meta$Type, Solute = s_meta$Solute)
ann_colors <- list(Type = c(Biosynthesis = "#440154FF", Transport = "#FDE725FF"),
                   Solute = c(Betaine = hue_pal()(length(levels(s_meta$Solute)))[1],
                              Cation = hue_pal()(length(levels(s_meta$Solute)))[2],
                              Ectoine = hue_pal()(length(levels(s_meta$Solute)))[3],
                              Glutamate = hue_pal()(length(levels(s_meta$Solute)))[4],
                              Glutamine = hue_pal()(length(levels(s_meta$Solute)))[5],
                              Hydroxyectoine = hue_pal()(length(levels(s_meta$Solute)))[6],
                              Proline = hue_pal()(length(levels(s_meta$Solute)))[7],
                              Trehalose = hue_pal()(length(levels(s_meta$Solute)))[8]))
phm2 <- pheatmap(s_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 fontsize = 10,
                 fontsize_row = 7,
                 cluster_rows = F,
                 cluster_cols = F,
                 gaps_row = c(9, 13, 18, 19, 21, 22, 26))
save_pheatmap_pdf(phm2, "SaltToleranceKO_genome.pdf")
# Note, this shows 34 KOs present in at least 1 genome. The list had 78, so 44 weren't present in any!



#### phnJ Genera ####
p <- read_xlsx("genomeSet86753_01-jul-2021.xlsx", sheet = 1) %>%
  group_by(Domain, Genome) %>%
  summarise(n = n()) %>%
  group_by(Domain) %>%
  summarise(n = n())

#### End Script ####