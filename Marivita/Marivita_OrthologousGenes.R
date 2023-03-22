# Orthologous protein analysis
# by Cliff Bueno de Mesquita, JGI, Spring 2021
# This script performs the analysis on orthologous proteins between 8 Marivita taxa
# The new focal taxa we are describing are referred to as MAG_28 and MAG_10
# The other 6 are sister taxa
# Takes the output from proteinortho as the input and uses ComplexUpset to plot



#### Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(ComplexUpset)
library(micropan)
library(magrittr)
`%notin%` <- Negate(`%in%`)

# Working directory
setwd("~/Desktop/Marivita/")



#### Proteinortho ####
# .tsv output file from proteinortho program
d <- read.delim("ProteinOrtho.proteinortho.tsv")

# Get sequences for ones not in the MAGs but in all others (use M. lacus)
# Then you can annotate these in Pfam to see what they are
others <- subset(d, X..Species == 5 & Marivita10_192.proteins.faa == "*" & Marivita28_82.proteins.faa == "*") %>%
  select(MarLac.proteins.faa) %>%
  separate(MarLac.proteins.faa, into = c("Seqs", "Repeats"), sep = ",") %>%
  select(-Repeats)
marlac <- readFasta("~/Desktop/Marivita/MarivitaProteins/MarLac.proteins.faa") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  select(-Junk) %>%
  filter(Header %in% others$Seqs)
# writeFasta(marlac, "OtherProts74.faa")

# Format data for upset plot
d[d == "*"] <- 0
d[d != 0] <- 1
d1 <- d %>%
  select(-X..Species, -Genes, -Alg..Conn.) %>%
  mutate(Num = seq(from = 1, to = nrow(d), by = 1),
         Protein = "Protein") %>%
  mutate(ID = paste(Protein, Num, sep = "")) %>%
  select(-Num, -Protein) %>%
  column_to_rownames(var = "ID") %>%
  mutate_if(is.character, as.integer)

# Make different data frame subsets
# In all others but MAG28
d2 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita28_82.proteins.faa == 0 & sum == 6) %>%
  select(-sum)

# In all others but MAG10
d3 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita10_192.proteins.faa == 0 & sum == 6) %>%
  select(-sum)

# In all others but MAG28 and MAG10
d4 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita28_82.proteins.faa == 0 & Marivita10_192.proteins.faa == 0 & sum == 5) %>%
  select(-sum)

# Only in MAG28 (none)
d5 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita28_82.proteins.faa == 1 & sum == 1) %>%
  select(-sum)

# Only in MAG10 (none)
d6 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita10_192.proteins.faa == 1 & sum == 1) %>%
  select(-sum)

# Only in both MAGs
d7 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(Marivita28_82.proteins.faa == 1 & Marivita10_192.proteins.faa == 1 & sum == 2) %>%
  select(-sum)

# All 8
d8 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 7) %>%
  select(-sum)

# Combine
d9 <- rbind(d2, d3, d4, d5, d6, d7, d8) %>%
  select(`Marivita10_192.proteins.faa`, `Marivita28_82.proteins.faa`, `MarCry.proteins.faa`,
         `MarGeo.proteins.faa`, `MarHal.proteins.faa`, `MarLZ.proteins.faa`, `MarLac.proteins.faa`)



#### KO ####
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

# In all others but MAG28
k2 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR1` == 0 & sum == 6) %>%
  select(-sum)

# In all others but MAG10
k3 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR2` == 0 & sum == 6) %>%
  select(-sum)

# In all others but MAG28 and MAG10
k4 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR1` == 0 & `M. sp. SBSPR2` == 0 & sum == 5) %>%
  select(-sum)

# Only in MAG28 (44)
k5 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR1` == 1 & sum == 1) %>%
  select(-sum)

# Only in MAG10 (92)
k6 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR2` == 1 & sum == 1) %>%
  select(-sum)

# Only in both MAGs (12)
k7 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(`M. sp. SBSPR1` == 1 & `M. sp. SBSPR2` == 1 & sum == 2) %>%
  select(-sum)

# All 8
k8 <- KO_table %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 7) %>%
  select(-sum)

# Combine
k9 <- rbind(k2, k3, k4, k5, k6, k7, k8) %>%
  select(`M. sp. SBSPR2`, `M. sp. SBSPR1`, `M. cryptomonadis`,
         `M. geojedonensis`, `M. hallyeonensis`, `M. sp. LZ-15-2`, `M. lacus`)
names(k9)[3] <- "M. cryptomonadis MP20-4"
names(k9)[6] <- "M. cryptomonadis LZ-15-2"
  
  
#### Complex Upset ####
# Single panel
# Shows: all species, all but focal species, pairwise with focal species, only in focal species
# Also add a labeling function to rename the genomes
# Also change the y axis to "Shared orthogroups", adjust margins, make names italic
n <- c("MarLac.proteins.faa" = "M. lacus",
       "MarLZ.proteins.faa" = "M. cryptomonadis LZ-15-2",
       "MarHal.proteins.faa" = "M. hallyeonensis",
       "MarGeo.proteins.faa" = "M. geojedonensis",
       "MarCry.proteins.faa" = "M. cryptomonadis MP20-4",
       "Marivita28_82.proteins.faa" = "M. sp. SBSPR1",
       "Marivita10_192.proteins.faa" = "M. sp. SBSPR2")

upset(d9,
      intersect = names(d9),
      set_sizes = FALSE,
      name = "Genomes",
      base_annotations = list('Shared orthogroups' = intersection_size(counts=T)),
      labeller = as_labeller(n),
      themes = upset_modify_themes(
        list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40),
                                                                       face = "italic"),
                                            plot.margin = margin(l = 30)))))

upset(k9,
      intersect = names(k9),
      set_sizes = FALSE,
      name = "Genomes",
      base_annotations = list('Shared KOs' = intersection_size(counts=T)),
      themes = upset_modify_themes(
        list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40),
                                                                       face = "italic"),
                                            plot.margin = margin(l = 30)))))



# Multipanel with proteinortho and KO results
# Don't show numbers - add them manually in powerpoint so they're bigger and all above bars
# (set counts = F)
pdf("~/Desktop/Marivita/Manuscript/UpsetMultipanel.pdf", width = 7, height = 5)
upset(d9,
      intersect = names(d9),
      wrap = TRUE,
      sort_sets = FALSE,
      set_sizes = FALSE,
      keep_empty_groups = TRUE,
      name = NULL,
      base_annotations = list('Shared orthogroups' = intersection_size(counts=F, text = aes(size = 2))),
      labeller = as_labeller(n),
      themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('a) proteinortho') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  upset(k9,
        intersect = names(k9),
        wrap = TRUE,
        sort_sets = FALSE,
        keep_empty_groups = TRUE,
        set_sizes = FALSE,
        name = NULL,
        base_annotations = list('Shared KOs' = intersection_size(counts=F, text = aes(size = 2))),
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('b) KEGG Orthology') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()



#### End script ####