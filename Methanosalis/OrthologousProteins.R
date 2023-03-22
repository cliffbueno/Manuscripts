# Orthologous protein analysis
# by Cliff Bueno de Mesquita, JGI, Spring 2021
# This script performs the analysis on orthologous proteins between 11 genomes
# The new focal taxon we are describing is referred to as MAG_48 (version with 48 contigs). 
# The others are 9 Methanolobus and 1 Methanomethylovorans
# Takes the output from proteinortho as the input and uses upset to plot
# Uses the ComplexUpset package, but the UpSetR package is another option

#### Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(ComplexUpset)
library(micropan)
`%notin%` <- Negate(`%in%`)

# Working directory
setwd("~/Desktop/Methanolobus/Proteinortho/")

# .tsv output file from proteinortho program
d <- read.delim("ProteinOrtho.proteinortho.tsv")

# Get sequences for ones not in MAG but in all others (use M. bombayensis)
# Then you can annotate these in Pfam to see what they are
others <- subset(d, X..Species == 10 & MAG_48.faa == "*") %>%
  select(Met_bom.faa) %>%
  separate(Met_bom.faa, into = c("Seqs", "Repeats"), sep = ",") %>%
  select(-Repeats)
metbom <- readFasta("~/Desktop/Methanolobus/MethanolobusProteins/Met_bom.faa") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  select(-Junk) %>%
  filter(Header %in% others$Seqs)
# writeFasta(metbom, "OtherProts85.faa") # Original
# writeFasta(metbom, "OtherProts61.faa") # Updated w Mm Hol

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
# Only containing MAG, pairs > 2
d2 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(MAG_48.faa == 1 & sum > 2) %>%
  select(-sum)

# In all others but MAG
d3 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(MAG_48.faa == 0 & sum == 10) %>%
  select(-sum)

# Only in MAG (0 by definition...)
d4 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(MAG_48.faa == 1 & sum == 1) %>%
  select(-sum)

# Pairwise with MAG
d5 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(MAG_48.faa == 1 & sum == 2) %>%
  select(-sum)

# All 11
d6 <- d1 %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 11) %>%
  select(-sum)

# Combine d6, d3, d5, d4 (in the example d4 is 0 and won't show anything)
d7 <- rbind(d6, d3, d5, d4) %>%
  select(Metv_hol.faa, MAG_48.faa, Met_psyp.faa, Met_zin.faa, Met_sy1.faa, Met_psyt.faa,
         Met_pro.faa, Met_bom.faa, Met_tin.faa, Met_vul1.faa, Met_vul.faa)



#### Complex Upset ####
# Default, all (too many)
upset(d1,
      intersect = names(d1),
      set_sizes = TRUE,
      name = "Genomes")

# Multipanel
upset(d2,
      intersect = names(d2),
      set_sizes = FALSE,
      name = "Genomes",
      min_size = 10,
      wrap = TRUE) +
  ggtitle('size > 10, sets > 2') +
  upset(d3,
        intersect = names(d3),
        set_sizes = FALSE,
        name = "Genomes",
        min_size = 10,
        wrap = TRUE) +
  ggtitle('Not in MAG') +
  upset(d5,
        intersect = names(d3),
        set_sizes = FALSE,
        name = "Genomes",
        wrap = TRUE) +
  ggtitle('Pairwise with MAG')

# Single panel
# Shows: all species, all but focal species, pairwise with focal species, only in focal species
# Also add a labeling function to rename the genomes
# Also change the y axis to "Shared orthogroups", adjust margins, make names italic
n <- c("Metv_hol.faa" = "Mmv. hollandica",
       "MAG_48.faa" = "Ms. sp. SBSPR1A",
       "Met_bom.faa" = "Ml. bombayensis",
       "Met_pro.faa" = "Ml. profundi",
       "Met_psyp.faa" = "Ml. psychrophilus",
       "Met_psyt.faa" = "Ml. psychrotolerans",
       "Met_sy1.faa" = "Ml. sp. SY-01",
       "Met_tin.faa" = "Ml. tindarius",
       "Met_vul.faa" = "Ml. vulcani B1d",
       "Met_vul1.faa" = "Ml. vulcani PL 12/M",
       "Met_zin.faa" = "Ml. zinderi")

upset(d7,
      intersect = names(d7),
      set_sizes = FALSE,
      name = "Genomes",
      base_annotations = list('Shared orthogroups' = intersection_size(counts=T)),
      labeller = as_labeller(n),
      themes = upset_modify_themes(
        list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40),
                                                                       face = "italic"),
                                            plot.margin = margin(l = 30)))))



# Multipanel with proteins and genes (use this) (need to make "k5" object in different script first)
pdf("~/Desktop/Methanolobus/Manuscript/Fig2.pdf", width = 7, height = 5)
upset(d7,
      intersect = names(d7),
      wrap = TRUE,
      sort_sets = FALSE,
      set_sizes = FALSE,
      keep_empty_groups = TRUE,
      name = NULL,
      base_annotations = list('Shared orthogroups' = intersection_size(counts=T, text = aes(size = 1.8))),
      labeller = as_labeller(n),
      themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('a) proteinortho') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  upset(k5,
        intersect = names(k5),
        wrap = TRUE,
        sort_sets = FALSE,
        keep_empty_groups = TRUE,
        set_sizes = FALSE,
        name = NULL,
        base_annotations = list('Shared KOs' = intersection_size(counts=T, text = aes(size = 1.8))),
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('b) KEGG Orthology') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()

# Revisions - reviewer wanted clearer text and numbers above
pdf("~/Desktop/Methanolobus/Manuscript/Upset_forPPT.pdf", width = 7, height = 5)
upset(d7,
      intersect = names(d7),
      height_ratio = 1,
      wrap = TRUE,
      sort_sets = FALSE,
      set_sizes = FALSE,
      keep_empty_groups = TRUE,
      name = NULL,
      base_annotations = list('Shared orthogroups' = intersection_size(counts=F, text = aes(size = 1.8))),
      labeller = as_labeller(n),
      themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), color = "black", size = 10), plot.margin = margin(l = 30))))) +
  ggtitle('a) proteinortho') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  upset(k5,
        intersect = names(k5),
        height_ratio = 1,
        wrap = TRUE,
        sort_sets = FALSE,
        keep_empty_groups = TRUE,
        set_sizes = FALSE,
        name = NULL,
        base_annotations = list('Shared KOs' = intersection_size(counts=F, text = aes(size = 1.8))),
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), color = "black", size = 10), plot.margin = margin(l = 30))))) +
  ggtitle('b) KEGG Orthology') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()

# Make png for google doc
png("~/Desktop/Methanolobus/Manuscript/Fig2.png", width = 7, height = 5, units = "in", res = 300)
upset(d7,
      intersect = names(d7),
      wrap = TRUE,
      sort_sets = FALSE,
      set_sizes = FALSE,
      keep_empty_groups = TRUE,
      name = NULL,
      base_annotations = list('Shared orthogroups' = intersection_size(counts=T, text = aes(size = 1.8))),
      labeller = as_labeller(n),
      themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('a) proteinortho') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  upset(k5,
        intersect = names(k5),
        wrap = TRUE,
        sort_sets = FALSE,
        keep_empty_groups = TRUE,
        set_sizes = FALSE,
        name = NULL,
        base_annotations = list('Shared KOs' = intersection_size(counts=T, text = aes(size = 1.8))),
        themes = upset_modify_themes(list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40), face = "italic"), plot.margin = margin(l = 30, b = -40))))) +
  ggtitle('b) KEGG Orthology') +
  theme(plot.title = element_text(hjust = 0.6, vjust = -0.5),
        plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()

#### End script ####