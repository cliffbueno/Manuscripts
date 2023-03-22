# Heatmaps for KOs involved in methanogenesis and salt tolerance
# Genome level - 11 genomes
# -MAG and 9 other Methanolobus genomes and 1 Methanomethylovorans genome
# 7 were on IMG, 4 others used BLAST Koala for KO assignment
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
show_col(hue_pal()(6))
show_col(viridis_pal()(3))
hue_pal()(13)[1]
viridis_pal()(3)[2]
save_pheatmap_pdf <- function(x, filename, width = 6.5, height = 6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
meta <- read_excel("~/Desktop/Metabolomics/Metabolites_10.20.20_Clean.xlsx", sheet = 1) %>%
  mutate(Site = as.factor(Site),
         Restoration = as.factor(Restoration))
setwd("~/Desktop/Methanolobus/Methanolobus_KO_profiles_IMG_KOALA/")

# KO table - import 6 from IMG and 4 from BLAST Koala
IMG_KO <- read.delim("stats_input", header = TRUE) %>%
  slice_tail(n = nrow(.)-1) %>%
  set_colnames(c("KO", "Ml. bombayensis", "Ml. psychrophilus", "Ml. profundi",
                 "Ml. tindarius", "Ml. zinderi", "Ml. vulcani PL 12/M")) %>%
  separate(KO, into = c("Junk", "KO"), sep = ":") %>%
  select(-Junk)

MAG_KO <- read.delim2("MAG_48_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MAG = n()) %>%
  set_colnames(c("KO", "Ms. sp. SBSPR1"))

MetPsyt_KO <- read.delim2("Met_psyt_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MetPsyt = n()) %>%
  set_colnames(c("KO", "Ml. psychrotolerans"))

MetSy1_KO <- read.delim2("Met_sy1_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MetSy1 = n()) %>%
  set_colnames(c("KO", "Ml. sp. SY-01"))

MetVul_KO <- read.delim2("Met_vul_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(MetVul = n()) %>%
  set_colnames(c("KO", "Ml. vulcani B1d"))

MetHol_KO <- read.delim2("methanomethylovorans_24-aug-2021/profile.txt", header = F) %>%
  filter(V1 == "Methanomethylovorans") %>%
  filter(V3 == 1) %>%
  select(-V1) %>%
  separate(V2, into = c("Junk", "KO"), sep = ":") %>%
  select(-Junk) %>%
  group_by(KO) %>%
  summarise(MetHol = n()) %>%
  set_colnames(c("KO", "Mmv. hollandica"))

# Merge into master table
KO_table <- full_join(IMG_KO, MAG_KO, by = "KO") %>%
  full_join(., MetPsyt_KO, by = "KO") %>%
  full_join(., MetSy1_KO, by = "KO") %>%
  full_join(., MetVul_KO, by = "KO") %>%
  full_join(., MetHol_KO, by = "KO") %>%
  column_to_rownames(var = "KO") %>%
  mutate_if(is.character, as.integer)

# Check copy number
to_check <- read_excel("~/Desktop/MetabolicModels/MethanogenesisRxn_for_Seed.xlsx", sheet = 3) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  filter(Pathway_General == "Methyl-CoM formation" |
           Pathway_General == "Methanogenesis")

KO_table_copy_check <- KO_table %>%
  rownames_to_column(var = "KO") %>%
  filter(KO %in% to_check$KEGG_KO)

# Make table of 1s and 0s
KO_table[is.na(KO_table)] <- 0
KO_table[KO_table > 0] <- 1
KO_table <- KO_table %>%
  rownames_to_column(var = "KO")


#### Methanogenesis ####
a <- read_excel("~/Desktop/MetabolicModels/KEGG_Modules.xlsx", sheet = 5) %>%
  filter(Type == "KO") %>%
  select(Type, KEGG, Formula, Module, Order) %>%
  mutate_if(is.character, as.factor) %>%
  group_by(KEGG) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Order)

meth_KOs <- KO_table %>%
  filter(KO %in% a$KEGG)

meth_meta <- left_join(meth_KOs, a, by = c("KO" = "KEGG")) %>%
  arrange(Order) %>%
  mutate(Module = factor(Module,
                         levels = c("Acetate", "Acetate-CO2", "CO2", "Methanol", "Methylamine", "All")))

# Optional - add KOs not present in KO table
# a$KEGG %notin% e$KO # 10 missing
# K00625, K00925, K00204, K13942, K14126, K14128, K22480, K22481, K22482, K22516 
# f <- as.data.frame(matrix(data = 0, nrow = 10, ncol = 7))
# names(f) <- names(e)
# f$KO <- c("K00625", "K00925", "K00204", "K13942", "K14126", 
#          "K14128", "K22480", "K22481", "K22482", "K22516")
# m <- rbind(e, f)

meth_KOs <- left_join(meth_KOs, a[,c(2,3,5)], by = c("KO" = "KEGG")) %>%
  mutate(KO_Def = paste(KO, Formula, sep = " ")) %>%
  arrange(Order) %>%
  select(-Order, -KO, -Formula) %>%
  select(`Ms. sp. SBSPR1`, `Ml. psychrophilus`, `Ml. zinderi`, `Ml. sp. SY-01`,
         `Ml. psychrotolerans`, `Ml. profundi`, `Ml. bombayensis`, `Ml. tindarius`,
         `Ml. vulcani PL 12/M`, `Ml. vulcani B1d` , KO_Def)

meth_mat <- meth_KOs %>%
  column_to_rownames(var = "KO_Def") %>%
  mutate_if(is.character, as.integer) %>%
  as.matrix()

# Pretty heatmap
setwd("~/Desktop/Methanolobus/Manuscript/")
ann_rows <- data.frame(row.names = rownames(meth_mat), Module = meth_meta$Module)
ann_colors <- list(Module = c(Acetate = "#F8766D", "Acetate-CO2" = "#B79F00", 
                              CO2 = "#00BA38", Methanol = "#00BFC4", 
                              Methylamine = "#619CFF", All = "#F564E3"))
phm1 <- pheatmap(meth_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 fontsize = 10,
                 fontsize_row = 6,
                 cluster_rows = F,
                 cluster_cols = F)
save_pheatmap_pdf(phm1, "MethaneModulesKO_genome.pdf")



#### _ Update ####
a <- read_excel("~/Desktop/MetabolicModels/MethanogenesisRxn_for_Seed.xlsx", sheet = 3) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order)

meth_KOs <- left_join(a, KO_table, by = c("KEGG_KO" = "KO")) %>%
  replace(is.na(.), 0) %>%
  select(8:19)
names(meth_KOs)[8] <- "Ms. sp. SBSPR1A"

meth_meta <- left_join(a, KO_table, by = c("KEGG_KO" = "KO")) %>%
  select(1:8)

meth_mat <- meth_KOs %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()

# Pretty heatmap
setwd("~/Desktop/Methanolobus/Manuscript/")
ann_rows <- data.frame(row.names = rownames(meth_mat), 
                       Pathway = meth_meta$Pathway_Specific,
                       Process = meth_meta$Pathway_General)
ann_colors <- list(Pathway = c(Acetate = brewer.pal(6, "Purples")[1], 
                               "Hydrogen/CO2" = brewer.pal(6, "Purples")[2], 
                               Trimethylamine = brewer.pal(6, "Purples")[3], 
                               Dimethylamine = brewer.pal(6, "Purples")[4], 
                               Methylamine = brewer.pal(6, "Purples")[5], 
                               Methanol = brewer.pal(6, "Purples")[6],
                               "Methyl-CoM reduction" = "#21908CFF",
                               Methanophenazine = brewer.pal(9, "YlOrRd")[1],
                               Ferredoxin = brewer.pal(9, "YlOrRd")[2],
                               "Ferredoxin-F420-H2" = brewer.pal(9, "YlOrRd")[3],
                               H2 = brewer.pal(9, "YlOrRd")[4]),
                   Process = c("Methyl-CoM formation" = "#440154FF",
                               "Methanogenesis" = "#21908CFF",
                               "CoB - CoM regeneration" = "#FDE725FF"))
phm1 <- pheatmap(meth_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 fontsize = 10,
                 fontsize_row = 5,
                 cluster_rows = F,
                 cluster_cols = F,
                 gaps_row = c(39, 42))
save_pheatmap_pdf(phm1, "Fig3.pdf")



#### _Precursors ####
b <- read_excel("~/Desktop/MetabolicModels/MethanogenesisRxn_for_Seed.xlsx", sheet = 4) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order)

prec_KOs <- left_join(b, KO_table, by = c("KEGG_KO" = "KO")) %>%
  replace(is.na(.), 0) %>%
  select(8:19)
names(prec_KOs)[8] <- "Ms. sp. SBSPR1A"

prec_meta <- left_join(b, KO_table, by = c("KEGG_KO" = "KO")) %>%
  select(1:8)

prec_mat <- prec_KOs %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()

# Pretty heatmap
setwd("~/Desktop/Methanolobus/Manuscript/")
ann_rows <- data.frame(row.names = rownames(prec_mat), 
                       Pathway = prec_meta$Pathway_Specific,
                       Process = prec_meta$Pathway_General)
ann_colors <- list(Pathway = c("Coenzyme B Synthesis" = "#440154FF", 
                               "Coenzyme M Synthesis I" = brewer.pal(9, "YlOrRd")[2],
                               "Coenzyme M Synthesis II" = brewer.pal(9, "YlOrRd")[4],
                               "Coenzyme M Synthesis" = brewer.pal(9, "YlOrRd")[6]),
                   Process = c("Coenzyme B" = "#440154FF",
                               "Coenzyme M" = "#FDE725FF"))
phm2 <- pheatmap(prec_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 fontsize = 10,
                 fontsize_row = 5,
                 cluster_rows = F,
                 cluster_cols = F,
                 gaps_row = c(9))
save_pheatmap_pdf(phm2, "FigS7.pdf")



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
  arrange(Order) %>%
  select(`Mmv. hollandica`, `Ms. sp. SBSPR1`, `Ml. psychrophilus`, `Ml. zinderi`, `Ml. sp. SY-01`,
         `Ml. psychrotolerans`, `Ml. profundi`, `Ml. bombayensis`, `Ml. tindarius`,
         `Ml. vulcani PL 12/M`, `Ml. vulcani B1d`, KO)
names(s)[2] <- "Ms. sp. SBSPR1A"
s_mat <- s %>%
  column_to_rownames(var = "KO") %>%
  mutate_if(is.character, as.integer) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(s_mat), Type = s_meta$Type, Solute = s_meta$Solute)
ann_colors <- list(Type = c(Biosynthesis = "#440154FF", Transport = "#FDE725FF"),
                   Solute = c(Betaine = hue_pal()(length(levels(s_meta$Solute)))[1],
                              Cation = hue_pal()(length(levels(s_meta$Solute)))[2],
                              Ectoine = hue_pal()(length(levels(s_meta$Solute)))[3],
                              Glutamate = hue_pal()(length(levels(s_meta$Solute)))[4],
                              Glutamine = hue_pal()(length(levels(s_meta$Solute)))[5],
                              Proline = hue_pal()(length(levels(s_meta$Solute)))[6],
                              Sucrose = hue_pal()(length(levels(s_meta$Solute)))[7]))
phm3 <- pheatmap(s_mat,
                 legend = F,
                 color = bluered(100),
                 border_color = NA,
                 scale = "none",
                 angle_col = 315,
                 annotation_row = ann_rows,
                 annotation_colors = ann_colors,
                 fontsize = 10,
                 fontsize_row = 6,
                 cluster_rows = F,
                 cluster_cols = F)
save_pheatmap_pdf(phm3, "Fig5.pdf")
# Note, this shows 22 KOs present in at least 1 genome. The list had 78, so 56 weren't present in any!



#### Orthologous Genes ####
# Make intersection plot for KOs
# Use object k5 in the OrthologousProteins.R script to make a multipanel graph with proteinrortho and KO results
library(ComplexUpset)

# Make different subsets
KO_table_pa <- KO_table %>%
  column_to_rownames(var = "KO")

# All 11
k1 <- KO_table_pa %>%
  mutate(sum = rowSums(.)) %>%
  filter(sum == 11) %>%
  select(-sum)

# In all 10 others but not MAG
k2 <- KO_table_pa %>%
  mutate(sum = rowSums(.)) %>%
  filter(`Ms. sp. SBSPR1` == 0 & sum == 10) %>%
  select(-sum)

# Pairwise with MAG
k3 <- KO_table_pa %>%
  mutate(sum = rowSums(.)) %>%
  filter(`Ms. sp. SBSPR1` == 1 & sum == 2) %>%
  select(-sum)

# Only in MAG (18)
k4 <- KO_table_pa %>%
  mutate(sum = rowSums(.)) %>%
  filter(`Ms. sp. SBSPR1` == 1 & sum == 1) %>%
  select(-sum)

k5 <- rbind(k1, k2, k3, k4) %>%
  select(`Mmv. hollandica`, `Ms. sp. SBSPR1`, `Ml. psychrophilus`, `Ml. zinderi`, 
         `Ml. sp. SY-01`, `Ml. psychrotolerans`, `Ml. profundi`, `Ml. bombayensis`,
         `Ml. tindarius`, `Ml. vulcani PL 12/M`, `Ml. vulcani B1d`)
names(k5)[2] <- "Ms. sp. SBSPR1A"

upset(k5,
      intersect = names(k5),
      set_sizes = FALSE,
      name = "Genomes",
      base_annotations = list('Shared KOs' = intersection_size(counts=T)),
      themes = upset_modify_themes(
        list('intersections_matrix' = theme(axis.text.y = element_text(margin = margin(l = -40),
                                                                       face = "italic"),
                                            plot.margin = margin(l = 30)))))

