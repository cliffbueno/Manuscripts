#### PanAf Chimpanzee Gut Microbiome Analysis ####
# Cleaned Analysis Script for PNAS Manuscript statistics and figures
# by Cliff Bueno de Mesquita, February-June 2020, cliff.buenodemesquita@colorado.edu
# Do all analyses on 560 samples with sex and no repeated sampling (except section 2S)

# Supplementary Figures:
# 1S. Alpha diversity rarefaction curves
# 2S. Within individual versus among individuals (within subspecies and site)
# 3S. Geography vs. Genetics, climate, vegetation
# 4S. Climate PCoA and Mantels
# 5S. Jaccard: Geography, Climate, Diet - all, elli, schwein, trog, verus
# 6S. Subspecies and diet (honey, termites) indicator taxa

# Main Figures:
# 1. Site Map (actually done in QGIS)
# 2. Taxa Drivers: a) Bacteria families by site/subsp b) Parasites by site/subsp
# 3. Multipanel PCoA (genetic, vegetation, 16s, 18s)
# 4. 16S - 18S Mantel tests and plot, all and each region.
# 5. BC: Geography, Vegetation - all, schwein, trog, verus
# 6. 16S BC & 18S Jac within vs. between sites a) elli, b) schwein, c) trog, d) verus
# 7. Diet (Algae, Honey, Nuts, Termites)

# Statistics
# 1. Permanovas
# 2. GDMs

# Note: subspecies names are replaced with West (verus), N-C (ellioti), Central (troglodytes), East (schweinfurthii) and multipanel figures will be ordered with regions going west to east left to right

######################################### Setup ###############################################
# Libraries
library(mctoolsr)
library(plyr)
library(tidyverse)
library(vegan)
library(RVAideMemoire)
library(reshape2)
library(cowplot)
library(ggtree)
library(BiodiversityR)
library(FSA)
library(TSdist)
library(philentropy)
library(lme4)
library(car)
library(phyloseq)
library(scales)
library(indicspecies)
library(fields)
library(geosphere)

# Functions
dist_geo <- function(lat_a, lon_a, lat_b, lon_b) { 
  if(anyNA(c(lat_a, lon_a, lat_b, lon_b))) return(NA) 
  round(distm(c(lon_a, lat_a), c(lon_b, lat_b), fun = distHaversine)/1000,2) 
}
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullp <- function(df) df[chull(df$Axis01p, df$Axis02p),]
find_hulld <- function(df) df[chull(df$Axis01d, df$Axis02d),]
find_hullc <- function(df) df[chull(df$Axis01c, df$Axis02c),]

# Environment
theme_set(theme_bw())
setwd("~/Desktop/CU/2Research/Chimp/Clean")
set.seed(500)

# Data
bac_tax_table_fp <- "../seqtab_wTax_mctoolsr_combine.txt"
bac_map_fp <- "../PanAF_Metadata_wClimData2.txt"
bac_input <- load_taxa_table(bac_tax_table_fp, bac_map_fp) # 1087
bac_input$map_loaded$sampleID <- rownames(bac_input$map_loaded)
euk_tax_table_fp <- "../18S/seqtab_wTax_mctoolsr.txt"
euk_map_fp <- "../PanAF_Metadata_wClimData2.txt"
euk_input <- load_taxa_table(euk_tax_table_fp, euk_map_fp) # 1120
euk_input$map_loaded$sampleID <- rownames(euk_input$map_loaded)
uniques <- read.csv("../uniques.csv") # 577 (430 w/1 sample, 147 w/> 1 samples & chose 1)
repeats <- read.csv("../repeats.csv") # 422
ad <- read.delim("observed_otus.txt") # 800
gen_mat_df <- read.delim("../AS Fst matrix 1_negatives_rmvd.txt", row.names = 1)
tree <- read.tree("../rep_abund_ubiq_phylo.tre") # Note tree unrooted and 3 archaea mixed in
abund_ubiq278 <- data.frame(asvID = tree$tip.label)
row.names(abund_ubiq278) <- abund_ubiq278$asvID
humans <- read.csv("../AmGut/human_geog_clim_bray.csv")
humans_sex <- read.csv("../AmGut/sex_human.csv")
humans_km <- readRDS("../AmGut/human_geogkm.rds")
plants <- read.csv("../habitat_indices_preliminary.csv")
plantspcomp <- read.csv("../species_compostion_matrix.csv", 
                        fileEncoding = "Latin1", check.names = F, row.names = 1)
biomAleman <- read.csv("../biom_Aleman.csv") %>%
  filter(Site_Name2 != "Delete") %>%
  select(Site_Name2, biome_Aleman)



################################## Data Processing ############################################
# Filtering and rarefying data
# Save objects as .rds
# To exactly reproduce the analysis, skip to next section and load .rds files
# Bacterial Dataset
# Remove blanks, controls and 'Do Not Use' samples
bac_input_filt <- filter_data(bac_input,
                         filter_cat = "Type", 
                         filter_vals = c("EXCLUDE","DoNotUse","Blank"))
# Remove contaminated samples
bac_input_filt <- filter_data(bac_input_filt,
                     filter_cat = "extr_name",
                     filter_vals = c("Gbo2-48", "Con1-19","Con2-27"))
# Remove mitochondrial, chloroplast and eukaryota sequences (376)
bac_input_filt <- filter_taxa_from_input(bac_input_filt, 
                                     taxa_to_remove=c("Chloroplast","Mitochondria", "Eukaryota")) 
# Remove reads that are unassigned at domain level (536)
bac_input_filt <- filter_taxa_from_input(bac_input_filt, at_spec_level = 2, 
                                         taxa_to_remove = "NA")

# Rarefy at 8000
bac_input_filt_rar <- single_rarefy(bac_input_filt, 8000) # 852 samples remaining
# Unique individuals with sex ID
bac_input_filt_rar_uniq <- filter_data(bac_input_filt_rar,
                              filter_cat = "sampleID",
                              keep_vals = uniques$sample_name)
bac_input_filt_rar_uniq <- filter_data(bac_input_filt_rar_uniq,
                                       filter_cat = "Sex",
                                       filter_vals = "*") # 560 samples
# saveRDS(bac_input_filt_rar_uniq, file = "bac_input_filt_rar_uniq.rds")
# Repeat individuals
bac_input_filt_rar_repeat <- filter_data(bac_input_filt_rar,
                                       filter_cat = "sampleID",
                                       keep_vals = repeats$sample_name) # 422 samples
# saveRDS(bac_input_filt_rar_repeat, file = "bac_input_filt_rar_repeat.rds")

### Eukaryote Dataset
# Remove blanks, controls and 'Do Not Use' samples
euk_input_filt <- filter_data(euk_input, 
                              filter_cat = "Type", 
                              filter_vals = c("EXCLUDE","DoNotUse","Blank","Water","NTC")) # 950
# Remove contaminated samples
euk_input_filt <- filter_data(euk_input_filt,
                              filter_cat = "extr_name",
                              filter_vals = c("Gbo2-48","Con1-19","Con2-27")) # 947
# Remove No Template Control NTC_1_b
euk_input_filt <- filter_data(euk_input_filt,
                              filter_cat = "sampleID",
                              filter_vals = "NTC_1_b") # 946
# Get rid of samples with less than 3100 seqs
euk_input_filt <- filter_data(euk_input_filt,
                              filter_cat = "sampleID",
                              filter_vals = c("Lib_y_013","TAI_R_410a",
                                              "Lib_C_433a","sob_399_a")) # 942
# Unique individuals with sex ID
euk_input_filt_uniq <- filter_data(euk_input_filt,
                                       filter_cat = "sampleID",
                                       keep_vals = uniques$sample_name)
euk_input_filt_uniq <- filter_data(euk_input_filt_uniq,
                                       filter_cat = "Sex",
                                       filter_vals = "*") # 560 samples
# Check number of samples with >= 50 reads
euk_input_filt_uniq_50 <- euk_input_filt_uniq
euk_input_filt_uniq_50$data_loaded$numsamps <- NA
for (i in 1:nrow(euk_input_filt_uniq_50$data_loaded)) {
  euk_input_filt_uniq_50$data_loaded$numsamps[i] <- sum(euk_input_filt_uniq_50$data_loaded[i,1:ncol(euk_input_filt_uniq_50$data_loaded)-1] > 49)
}
head(euk_input_filt_uniq_50$data_loaded$numsamps)
euk_input_filt_uniq_50$taxonomy_loaded$numsamps <- euk_input_filt_uniq_50$data_loaded$numsamps
View(euk_input_filt_uniq_50$taxonomy_loaded)

# Keep ubiquitous known parasites
# The following 11 ASV's present with >= 50 reads in 10% of samples, plus Blastocystis, Strongyloides, and 2 Entamoeba present with >= 50 reads in >5% of samples.
use <- c("ASV_1","ASV_11","ASV_12","ASV_14","ASV_15","ASV_2","ASV_22","ASV_26","ASV_31","ASV_51","ASV_8","ASV_59","ASV_90","ASV_80","ASV_70")
euk_input_filt_uniq_tenpercbs <- filter_taxa_from_input(euk_input_filt_uniq,
                                                        taxa_IDs_to_keep = use)
# Update taxonomy
# write.csv(euk_input_filt_uniq_tenpercbs$taxonomy_loaded, "parasite_taxonomy.csv")
# For this, stay conservative. Then in table can add some BLAST results
# View(euk_input_filt_uniq_tenpercbs$taxonomy_loaded
euk_input_filt_uniq_tenpercbs$taxonomy_loaded$taxonomy10 <- c("Troglodytella_abbrassarti",
                                                              "Dientamoeba_fragilis",
                                                              "Trichomonadidae_12",
                                                              "Haemonchus_contortus",
                                                              "Trichomonadidae_15",
                                                              "Tetratrichomonas_2",
                                                              "Blepharocorys_uncinata",
                                                              "Chilomastix_mesnili",
                                                              "Tetratrichomonas_31",
                                                              "Trichomonadidae_51",
                                                              "Blastocystis_59",
                                                              "Entamoeba_muris",
                                                              "Trichomonadidae_8",
                                                              "Entamoeba_hartmanni",
                                                              "Strongyloides_fuelleborni")
### Combined datasets
# New list
bac_euk <- bac_input_filt_rar_uniq
# Make month factor
bac_euk$map_loaded$collectionmonth <- as.factor(bac_euk$map_loaded$collectionmonth)
# Diet Column
bac_euk$map_loaded$Diet <- paste0(bac_euk$map_loaded$algae,
                                  bac_euk$map_loaded$ants,
                                  bac_euk$map_loaded$fruit,
                                  bac_euk$map_loaded$honey,
                                  bac_euk$map_loaded$marrow,
                                  bac_euk$map_loaded$meat,
                                  bac_euk$map_loaded$nuts,
                                  bac_euk$map_loaded$`palm heart`,
                                  bac_euk$map_loaded$termites,
                                  bac_euk$map_loaded$tubers,
                                  bac_euk$map_loaded$water)
bac_euk$map_loaded$Diet <- as.factor(bac_euk$map_loaded$Diet)
# Update site names
levels(bac_euk$map_loaded$Site)
bac_euk$map_loaded$Site <- revalue(bac_euk$map_loaded$Site, 
                                   c("Ugalla"="Issa",
                                     "Comoe_Geprenaf"="Comoe_GEPRENAF",
                                     "Tai_E"="Tai_Ecotourism",
                                     "Tai_R"="Tai_Recherche"))
levels(bac_euk$map_loaded$Site)
# Reorder sites for graphs (sorted by subspecies, country, alphabet)
bac_euk$map_loaded$Site <- factor(bac_euk$map_loaded$Site, 
                                  levels = c("Mt_Cameroon","Gashaka",
                                             "Gishwati","Nyungwe","Issa","Budongo","Bwindi",
                                             "Conkouati","Goualougo","Loango","Lope","Mts_de_Cristal",
                                             "Bakoun","Sangaredi","Sobeya","Sobory","Boe","Comoe_GEPRENAF","Djouroutou","Mt_Sangbe","Tai_Ecotourism","Tai_Recherche","East_Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi"))
# Add eukaryotes. Very important to line up dataframes!!
bac_euk$map_loaded$sampleID == euk_input_filt_uniq_tenpercbs$map_loaded$sampleID
rownames(bac_euk$map_loaded) == colnames(bac_euk$data_loaded)
rownames(bac_euk$data_loaded) == rownames(bac_euk$taxonomy_loaded)
bac_euk$map_loaded <- bac_euk$map_loaded[order(bac_euk$map_loaded$sampleID),]
rownames(bac_euk$map_loaded) == colnames(bac_euk$data_loaded) # Bad! 
bac_euk$data_loaded <- bac_euk$data_loaded[,order(names(bac_euk$data_loaded))]
rownames(bac_euk$map_loaded) == colnames(bac_euk$data_loaded) # Good!
rownames(bac_euk$data_loaded) == rownames(bac_euk$taxonomy_loaded) # Good!
euk_input_filt_uniq_tenpercbs$map_loaded <- euk_input_filt_uniq_tenpercbs$map_loaded[order(euk_input_filt_uniq_tenpercbs$map_loaded$sampleID),]
euk_input_filt_uniq_tenpercbs$data_loaded <- euk_input_filt_uniq_tenpercbs$data_loaded[,order(names(euk_input_filt_uniq_tenpercbs$data_loaded))]
rownames(euk_input_filt_uniq_tenpercbs$map_loaded) == colnames(euk_input_filt_uniq_tenpercbs$data_loaded) # Good!
rownames(euk_input_filt_uniq_tenpercbs$data_loaded) == rownames(euk_input_filt_uniq_tenpercbs$taxonomy_loaded) # Good!
bac_euk$map_loaded$sampleID == euk_input_filt_uniq_tenpercbs$map_loaded$sampleID # Good!
# Add euks. Combine ASV 12 and 51 because highly correlated r = 0.82
# Store as separate dataframes to add to bacterial metadata
ASV_8 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                taxa_IDs_to_keep = "ASV_8")
ASV_12 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_12")
ASV_15 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_15")
ASV_51 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_51")
ASV_2 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                taxa_IDs_to_keep = "ASV_2")
ASV_31 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_31")
ASV_11 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_11")
ASV_14 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_14")
ASV_59 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_59")
ASV_90 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_90")
ASV_1 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                taxa_IDs_to_keep = "ASV_1")
ASV_22 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_22")
ASV_26 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_26")
ASV_70 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_70")
ASV_80 <- filter_taxa_from_input(euk_input_filt_uniq_tenpercbs,
                                 taxa_IDs_to_keep = "ASV_80")
bac_euk$map_loaded$sampleID == ASV_1$map_loaded$sampleID
bac_euk$map_loaded$Blastocystis_59 <- as.vector(t(ASV_59$data_loaded))
bac_euk$map_loaded$Blepharocorys_uncinata <- as.vector(t(ASV_22$data_loaded))
bac_euk$map_loaded$Chilomastix_mesnili <- as.vector(t(ASV_26$data_loaded))
bac_euk$map_loaded$Dientamoeba_fragilis <- as.vector(t(ASV_11$data_loaded))
bac_euk$map_loaded$Entamoeba_hartmanni <- as.vector(t(ASV_80$data_loaded))
bac_euk$map_loaded$Entamoeba_muris <- as.vector(t(ASV_70$data_loaded))
bac_euk$map_loaded$Haemonchus_contortus <- as.vector(t(ASV_14$data_loaded))
bac_euk$map_loaded$Strongyloides_fuelleborni <- as.vector(t(ASV_90$data_loaded))
bac_euk$map_loaded$Tetratrichomonas_2 <- as.vector(t(ASV_2$data_loaded))
bac_euk$map_loaded$Tetratrichomonas_31 <- as.vector(t(ASV_31$data_loaded))
bac_euk$map_loaded$Trichomonadidae_12_51<-as.vector(t(ASV_12$data_loaded)+t(ASV_51$data_loaded))
bac_euk$map_loaded$Trichomonadidae_15 <- as.vector(t(ASV_15$data_loaded))
bac_euk$map_loaded$Trichomonadidae_8 <- as.vector(t(ASV_8$data_loaded))
bac_euk$map_loaded$Troglodytella_abrassarti <- as.vector(t(ASV_1$data_loaded))
# Parasite Presence/Absence (>= 50 to 1, < 50 to 0)
bac_euk$map_loaded <- mutate(bac_euk$map_loaded, 
                             Blastocystis_59 = ifelse(Blastocystis_59 >= 50, 1, 0),
                             Blepharocorys_uncinata = ifelse(Blepharocorys_uncinata >= 50, 1, 0),
                             Chilomastix_mesnili = ifelse(Chilomastix_mesnili >= 50, 1, 0),
                             Dientamoeba_fragilis = ifelse(Dientamoeba_fragilis >= 50, 1, 0),
                             Entamoeba_hartmanni = ifelse(Entamoeba_hartmanni >= 50, 1, 0),
                             Entamoeba_muris = ifelse(Entamoeba_muris >= 50, 1, 0),
                             Haemonchus_contortus = ifelse(Haemonchus_contortus >= 50, 1, 0),
                             Strongyloides_fuelleborni=ifelse(Strongyloides_fuelleborni>=50,1,0),
                             Tetratrichomonas_2 = ifelse(Tetratrichomonas_2 >= 50, 1, 0),
                             Tetratrichomonas_31 = ifelse(Tetratrichomonas_31 >= 50, 1, 0),
                             Trichomonadidae_12_51 = ifelse(Trichomonadidae_12_51 >= 50, 1, 0),
                             Trichomonadidae_15 = ifelse(Trichomonadidae_15 >= 50, 1, 0),
                             Trichomonadidae_8 = ifelse(Trichomonadidae_8 >= 50, 1, 0),
                             Troglodytella_abrassarti = ifelse(Troglodytella_abrassarti>=50,1,0))

# Diet Presence Absence (y to 1, n to 0)
bac_euk$map_loaded <- mutate(bac_euk$map_loaded,
                             algae = ifelse(algae == "y", 1, 0),
                             ants = ifelse(ants == "y", 1, 0),
                             fruit = ifelse(fruit == "y", 1, 0),
                             honey = ifelse(honey == "y", 1, 0),
                             marrow = ifelse(marrow == "y", 1, 0),
                             meat = ifelse(meat == "y", 1, 0),
                             nuts = ifelse(nuts == "y", 1, 0),
                             `palm heart` = ifelse(`palm heart` == "y", 1, 0),
                             termites = ifelse(termites == "y", 1, 0),
                             tubers = ifelse(tubers == "y", 1, 0),
                             water = ifelse(water == "y", 1, 0))
# This removed rownames so readd
rownames(bac_euk$map_loaded) <- bac_euk$map_loaded$sampleID
rownames(bac_euk$map_loaded) == colnames(bac_euk$data_loaded)

# Make protein column for any diet containing ants, termites or marrow
bac_euk$map_loaded$Protein <- NA
for (i in 1:nrow(bac_euk$map_loaded)) {
  if (bac_euk$map_loaded$ants[i] == 1 |
      bac_euk$map_loaded$termites[i] == 1 |
      bac_euk$map_loaded$marrow[i] == 1) {
    bac_euk$map_loaded$Protein[i] <- "High"
  } else {bac_euk$map_loaded$Protein[i] = "Low"}
}
class(bac_euk$map_loaded$Protein)
bac_euk$map_loaded$Protein <- as.factor(bac_euk$map_loaded$Protein)
table(bac_euk$map_loaded$Protein)
table(bac_euk$map_loaded$Sex)
# saveRDS(bac_euk, file = "bac_euk.rds")



#################################### Start Here ###############################################
bac_input_filt_rar_uniq <- readRDS("bac_input_filt_rar_uniq.rds")
bac_input_filt_rar_repeat <- readRDS("bac_input_filt_rar_repeat.rds")
bac_euk <- readRDS("bac_euk.rds")
# Sample ID list of the 560
samp_list <- bac_euk$map_loaded %>%
  select(sampleID)
sum(duplicated(samp_list$sampleID))
# write.csv(samp_list, "microbiome_samples_560.csv", row.names = F)
rn <- rownames(bac_euk$map_loaded)
bac_euk$map_loaded <- left_join(bac_euk$map_loaded, plants, by = c("Site" = "sites"))
rownames(bac_euk$map_loaded) <- rn
# Change Tai honey and termites to 1
for (i in 1:nrow(bac_euk$map_loaded)) {
  if (bac_euk$map_loaded$Site[i] == "Tai_Ecotourism" |
      bac_euk$map_loaded$Site[i] == "Tai_Recherche") {
    bac_euk$map_loaded$honey[i] <- 1
  }
}
for (i in 1:nrow(bac_euk$map_loaded)) {
  if (bac_euk$map_loaded$Site[i] == "Tai_Ecotourism" |
      bac_euk$map_loaded$Site[i] == "Tai_Recherche") {
    bac_euk$map_loaded$termites[i] <- 1
  }
}

# Subset by subspecies
elli <- filter_data(bac_euk,
                    filter_cat = "subspecies",
                    keep_vals = "ellioti") # 28
table(elli$map_loaded$Site)
table(elli$map_loaded$Protein)
table(elli$map_loaded$collectionmonth)
schwein <- filter_data(bac_euk,
                       filter_cat = "subspecies",
                       keep_vals = "schweinfurthii") # 134
table(schwein$map_loaded$Site)
table(schwein$map_loaded$Protein)
table(schwein$map_loaded$collectionmonth)
schwein$map_loaded$Diet <- revalue(schwein$map_loaded$Diet, 
                                   c("nnnnnnnnnnn" = "Nothing",
                                     "nnnynnnnnnn" = "Honey",
                                     "nynnnnnnyyn" = "Ants, Termites, Tubers",
                                     "nynynnnnnnn" = "Honey, Ants"))
trog <- filter_data(bac_euk, 
                    filter_cat = "subspecies",
                    keep_vals = "troglodytes") # 86
table(trog$map_loaded$Site)
table(trog$map_loaded$Protein)
table(trog$map_loaded$collectionmonth)
trog$map_loaded$Diet <- revalue(trog$map_loaded$Diet, 
                                c("nnnnnnnnynn" = "Termites",
                                  "nnnynnnnnnn" = "Honey",
                                  "nynynnnnnny" = "Ants, Honey, Water",
                                  "nynynnnnyny" = "Ants, Honey, Water, Termites",
                                  "nynyynnnyny" = "Ants, Honey, Water, Termites, Marrow"))
verus <- filter_data(bac_euk, 
                     filter_cat = "subspecies",
                     keep_vals = "verus") # 312
table(verus$map_loaded$Site)
table(verus$map_loaded$Protein)
table(verus$map_loaded$collectionmonth)

# Get site table with countries, climate, tree info
site_info <- select(bac_euk$map_loaded,
                    subspecies, Country, biom_AK, Site, AnTemp, Seasnlty, AnPercip, PrcpSeasnlty,
                    habitat, tree_density_ha, basal_area_m2_ha, species_richness,
                    Human_Footprint, lat, long)
site_info <- site_info[ !duplicated(site_info$Site), ]
# write.csv(site_info, "site_info.csv", row.names = FALSE)



### Matrices, all sample level
# Geography Euclidean (changed to km)
# For degrees use geography.distance <- dist(geography, method = "euclidean")
geography <- select(bac_euk$map_loaded, 
                    long, lat)
geography.distance.mat <- rdist.earth(geography, miles = FALSE, R = 6371)
rownames(geography.distance.mat) <- colnames(geography.distance.mat) <- rownames(geography.distance.mat)
geography.distance <- as.dist(geography.distance.mat)
elli.geography <- select(elli$map_loaded, 
                         long, lat)
elli.geography.distance.mat <- rdist.earth(elli.geography, miles = FALSE, R = 6371)
rownames(elli.geography.distance.mat) <- colnames(elli.geography.distance.mat) <- rownames(elli.geography.distance.mat)
elli.geography.distance <- as.dist(elli.geography.distance.mat)
schwein.geography <- select(schwein$map_loaded, 
                            long, lat)
schwein.geography.distance.mat <- rdist.earth(schwein.geography, miles = FALSE, R = 6371)
rownames(schwein.geography.distance.mat) <- colnames(schwein.geography.distance.mat) <- rownames(schwein.geography.distance.mat)
schwein.geography.distance <- as.dist(schwein.geography.distance.mat)
trog.geography <- select(trog$map_loaded, 
                         long, lat)
trog.geography.distance.mat <- rdist.earth(trog.geography, miles = FALSE, R = 6371)
rownames(trog.geography.distance.mat) <- colnames(trog.geography.distance.mat) <- rownames(trog.geography.distance.mat)
trog.geography.distance <- as.dist(trog.geography.distance.mat)
verus.geography <- select(verus$map_loaded, 
                          long, lat)
verus.geography.distance.mat <- rdist.earth(verus.geography, miles = FALSE, R = 6371)
rownames(verus.geography.distance.mat) <- colnames(verus.geography.distance.mat) <- rownames(verus.geography.distance.mat)
verus.geography.distance <- as.dist(verus.geography.distance.mat)

# Climate Euclidean (scale variables because not same unit!)
climate <- select(bac_euk$map_loaded, 
                  AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
climate <- as.data.frame(scale(climate))
climate.distance <- dist(climate, method = "euclidean")
elli.climate <- select(elli$map_loaded, 
                       AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
elli.climate <- as.data.frame(scale(elli.climate))
elli.climate.distance <- dist(elli.climate, method = "euclidean")
schwein.climate <- select(schwein$map_loaded,
                          AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
schwein.climate <- as.data.frame(scale(schwein.climate))
schwein.climate.distance <- dist(schwein.climate, method = "euclidean")
trog.climate <- select(trog$map_loaded, 
                       AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
trog.climate <- as.data.frame(scale(trog.climate))
trog.climate.distance <- dist(trog.climate, method = "euclidean")
verus.climate <- select(verus$map_loaded, 
                        AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
verus.climate <- as.data.frame(scale(verus.climate))
verus.climate.distance <- dist(verus.climate, method = "euclidean")

# Diet Jaccard. If warning given, must change NA to 0 or 1 (for sites with all zeroes)
diet <- select(bac_euk$map_loaded, 
               algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
               termites, tubers, water)
diet.distance <- vegdist(diet, method = "jaccard")
diet.distance <- dist.zeroes(diet, diet.distance)
elli.diet <- select(elli$map_loaded, 
                    algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                    termites, tubers, water)
elli.diet.distance <- vegdist(elli.diet, method = "jaccard")
schwein.diet <- select(schwein$map_loaded,
                       algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                       termites, tubers, water)
schwein.diet.distance <- vegdist(schwein.diet, method = "jaccard")
schwein.diet.distance <- dist.zeroes(schwein.diet, schwein.diet.distance)
trog.diet <- select(trog$map_loaded, 
                    algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                    termites, tubers, water)
trog.diet.distance <- vegdist(trog.diet, method = "jaccard")
verus.diet <- select(verus$map_loaded, 
                     algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                     termites, tubers, water)
verus.diet.distance <- vegdist(verus.diet, method = "jaccard")
verus.diet.distance <- dist.zeroes(verus.diet, verus.diet.distance)

# Parasite Jaccard. If warning given, must change NA to 0 or 1 (for sites with all zeroes)
parasites <- select(bac_euk$map_loaded,
                    Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                    Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                    Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                    Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                    Trichomonadidae_8, Troglodytella_abrassarti)
parasites.distance <- vegdist(parasites, method = "jaccard")
parasites.distance <- dist.zeroes(parasites, parasites.distance)
elli.parasites <- select(elli$map_loaded, 
                         Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                         Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                         Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                         Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                         Trichomonadidae_8, Troglodytella_abrassarti)
elli.parasites.distance <- vegdist(elli.parasites, method = "jaccard")
elli.parasites.distance <- dist.zeroes(elli.parasites, elli.parasites.distance)
schwein.parasites <- select(schwein$map_loaded,
                            Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                            Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                            Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                            Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                            Trichomonadidae_8, Troglodytella_abrassarti)
schwein.parasites.distance <- vegdist(schwein.parasites, method = "jaccard")
trog.parasites <- select(trog$map_loaded, 
                         Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                         Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                         Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                         Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                         Trichomonadidae_8, Troglodytella_abrassarti)
trog.parasites.distance <- vegdist(trog.parasites, method = "jaccard")
verus.parasites <- select(verus$map_loaded, 
                          Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                          Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                          Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                          Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                          Trichomonadidae_8, Troglodytella_abrassarti)
verus.parasites.distance <- vegdist(verus.parasites, method = "jaccard")

# Bacteria Bray
bacteria.distance <- calc_dm(bac_euk$data_loaded, method = "bray_sq_trans")
elli.bacteria.distance <- calc_dm(elli$data_loaded, method = "bray_sq_trans")
schwein.bacteria.distance <- calc_dm(schwein$data_loaded, method = "bray_sq_trans")
trog.bacteria.distance <- calc_dm(trog$data_loaded, method = "bray_sq_trans")
verus.bacteria.distance <- calc_dm(verus$data_loaded, method = "bray_sq_trans")
bac_hel <- decostand(t(bac_euk$data_loaded), method = "hellinger")



### Site level matrices (Geography, Climate, Diet, Genetics). These are site level
# Geography Euclidean
sl.geography <- select(bac_euk$map_loaded,
                       long, lat)
sl.geography <- sl.geography[ !duplicated(sl.geography$long), ]
sl.geography.distance.mat <- rdist.earth(sl.geography, miles = FALSE, R = 6371)
rownames(sl.geography.distance.mat) <- colnames(sl.geography.distance.mat) <- rownames(sl.geography.distance.mat)
sl.geography.distance <- as.dist(sl.geography.distance.mat)
sl.elli.geography <- select(elli$map_loaded, 
                            long, lat)
sl.elli.geography <- sl.elli.geography[ !duplicated(sl.elli.geography$long), ]
sl.elli.geography.distance.mat <- rdist.earth(sl.elli.geography, miles = FALSE, R = 6371)
rownames(sl.elli.geography.distance.mat) <- colnames(sl.elli.geography.distance.mat) <- rownames(sl.elli.geography.distance.mat)
sl.elli.geography.distance <- as.dist(sl.elli.geography.distance.mat)
sl.schwein.geography <- select(schwein$map_loaded, 
                               long, lat)
sl.schwein.geography <- sl.schwein.geography[ !duplicated(sl.schwein.geography$long), ]
sl.schwein.geography.distance.mat <- rdist.earth(sl.schwein.geography, miles = FALSE, R = 6371)
rownames(sl.schwein.geography.distance.mat) <- colnames(sl.schwein.geography.distance.mat) <- rownames(sl.schwein.geography.distance.mat)
sl.schwein.geography.distance <- as.dist(sl.schwein.geography.distance.mat)
sl.trog.geography <- select(trog$map_loaded, 
                            long, lat)
sl.trog.geography <- sl.trog.geography[ !duplicated(sl.trog.geography$long), ]
sl.trog.geography.distance.mat <- rdist.earth(sl.trog.geography, miles = FALSE, R = 6371)
rownames(sl.trog.geography.distance.mat) <- colnames(sl.trog.geography.distance.mat) <- rownames(sl.trog.geography.distance.mat)
sl.trog.geography.distance <- as.dist(sl.trog.geography.distance.mat)
sl.verus.geography <- select(verus$map_loaded, 
                             long, lat)
sl.verus.geography <- sl.verus.geography[ !duplicated(sl.verus.geography$long), ]
sl.verus.geography.distance.mat <- rdist.earth(sl.verus.geography, miles = FALSE, R = 6371)
rownames(sl.verus.geography.distance.mat) <- colnames(sl.verus.geography.distance.mat) <- rownames(sl.verus.geography.distance.mat)
sl.verus.geography.distance <- as.dist(sl.verus.geography.distance.mat)

# Climate Euclidean (scale variables because not same unit!). Collapse by Site.
sl.climate <- select(bac_euk$map_loaded, 
                     AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site)
sl.climate <- sl.climate[ !duplicated(sl.climate$Site), ]
sl.climate <- select(sl.climate, -Site)
sl.climate <- as.data.frame(scale(sl.climate))
sl.climate.distance <- dist(sl.climate, method = "euclidean")
sl.elli.climate <- select(elli$map_loaded, 
                          AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site)
sl.elli.climate <- sl.elli.climate[ !duplicated(sl.elli.climate$Site), ]
sl.elli.climate <- select(sl.elli.climate, -Site)
sl.elli.climate <- as.data.frame(scale(sl.elli.climate))
sl.elli.climate.distance <- dist(sl.elli.climate, method = "euclidean")
sl.schwein.climate <- select(schwein$map_loaded,
                          AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site)
sl.schwein.climate <- sl.schwein.climate[ !duplicated(sl.schwein.climate$Site), ]
sl.schwein.climate <- select(sl.schwein.climate, -Site)
sl.schwein.climate <- as.data.frame(scale(sl.schwein.climate))
sl.schwein.climate.distance <- dist(sl.schwein.climate, method = "euclidean")
sl.trog.climate <- select(trog$map_loaded, 
                          AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site)
sl.trog.climate <- sl.trog.climate[ !duplicated(sl.trog.climate$Site), ]
sl.trog.climate <- select(sl.trog.climate, -Site)
sl.trog.climate <- as.data.frame(scale(sl.trog.climate))
sl.trog.climate.distance <- dist(sl.trog.climate, method = "euclidean")
sl.verus.climate <- select(verus$map_loaded, 
                        AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site)
sl.verus.climate <- sl.verus.climate[ !duplicated(sl.verus.climate$Site), ]
sl.verus.climate <- select(sl.verus.climate, -Site)
sl.verus.climate <- as.data.frame(scale(sl.verus.climate))
sl.verus.climate.distance <- dist(sl.verus.climate, method = "euclidean")

# Diet Jaccard. If warning given, must change NA to 0 or 1 (for sites with all zeroes)
sl.diet <- select(bac_euk$map_loaded, 
                  algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                  termites, tubers, water, Site)
sl.diet <- sl.diet[ !duplicated(sl.diet$Site), ]
sl.diet <- select(sl.diet, -Site)
sl.diet.distance <- vegdist(sl.diet, method = "jaccard")
sl.diet.distance <- dist.zeroes(sl.diet, sl.diet.distance)
sl.elli.diet <- select(elli$map_loaded, 
                       algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                       termites, tubers, water, Site)
sl.elli.diet <- sl.elli.diet[ !duplicated(sl.elli.diet$Site), ]
sl.elli.diet <- select(sl.elli.diet, -Site)
sl.elli.diet.distance <- vegdist(sl.elli.diet, method = "jaccard")
sl.schwein.diet <- select(schwein$map_loaded,
                          algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                          termites, tubers, water, Site)
sl.schwein.diet <- sl.schwein.diet[ !duplicated(sl.schwein.diet$Site), ]
sl.schwein.diet <- select(sl.schwein.diet, -Site)
sl.schwein.diet.distance <- vegdist(sl.schwein.diet, method = "jaccard")
sl.trog.diet <- select(trog$map_loaded, 
                       algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                       termites, tubers, water, Site)
sl.trog.diet <- sl.trog.diet[ !duplicated(sl.trog.diet$Site), ]
sl.trog.diet <- select(sl.trog.diet, -Site)
sl.trog.diet.distance <- vegdist(sl.trog.diet, method = "jaccard")
sl.verus.diet <- select(verus$map_loaded, 
                        algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                        termites, tubers, water, Site)
sl.verus.diet <- sl.verus.diet[ !duplicated(sl.verus.diet$Site), ]
sl.verus.diet <- select(sl.verus.diet, -Site)
sl.verus.diet.distance <- vegdist(sl.verus.diet, method = "jaccard")

# Write file for the 560 filtered samples data, for alpha rarefaction
# write.csv(bac_euk$data_loaded, "seqtab560samplesrarefied8000.csv")

# Habitat/Vegetation Data
# Habitat categories
# outlier is Outamba
ggplot(data = site_info, aes(habitat, basal_area_m2_ha)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.5) +
  labs(x = "Biome",
       y = "Basal Area (m2/ha)") +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))

# Climate distance vs. Plant Bray-Curtis
# Climate distance, no Goualougo or Lope, match with plant matrix.
sl.climate.h <- select(bac_euk$map_loaded, 
                       AnTemp, AnPercip, Seasnlty, PrcpSeasnlty, Site, subspecies)
sl.climate.h <- sl.climate.h[ !duplicated(sl.climate.h$Site), ]
sl.climate.h <- subset(sl.climate.h, Site != "Goualougo")
sl.climate.h <- subset(sl.climate.h, Site != "Lope")
sl.climate.h <- sl.climate.h[order(sl.climate.h$Site),]
# Check match
rownames(plantspcomp) == sl.climate.h$Site
plantmeta <- sl.climate.h
sl.climate.h <- select(sl.climate.h, -Site, -subspecies)
sl.climate.h <- as.data.frame(scale(sl.climate.h))
sl.climate.distance.h <- dist(sl.climate.h, method = "euclidean")
# Plant distance
plantrel <- decostand(plantspcomp, method = "total")
plantbray.dist <- vegdist(plantrel, method = "bray")
mantel(sl.climate.distance.h, plantbray.dist) # r = 0.54, p = 0.001
sl.climate.distance.h <- as.matrix(sl.climate.distance.h)
sl.climate.distance.h[upper.tri(sl.climate.distance.h, diag = TRUE)] <- NA
plantbray <- as.matrix(plantbray.dist)
plantbray[upper.tri(plantbray, diag = TRUE)] <- NA
plant_clim_dist <- data.frame('clim' = na.omit(c(sl.climate.distance.h)),
                              'plants' = na.omit(c(plantbray)))
label.df.pcd <- data.frame(x = 3,
                           y = 0.125,
                           label = c("n = 27, r = 0.54, p = 0.001"))
pdf(file = "Veg_Clim.pdf", width = 6.5, height = 4)
ggplot(plant_clim_dist, aes(clim, plants)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  labs(x = "Climate Distance (Euclidean)",
       y = "Plant Dissimilarity (Bray-Curtis)") +
  geom_text(data = label.df.pcd, aes(x = x, y = y, label = label), size = 4) +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
dev.off()

# Remove Goualougo and Lope, add veg data, calc bacteria and plant dissimilary matrices
bac_euk_veg <- filter_data(bac_euk, 
                           filter_cat = "Site", 
                           filter_vals = c("Goualougo","Lope"))
plantspcompjoin <- plantspcomp
plantspcompjoin$Site <- rownames(plantspcompjoin)
rn3 <- rownames(bac_euk_veg$map_loaded)
bac_euk_veg$map_loaded <- left_join(bac_euk_veg$map_loaded, plantspcompjoin, 
                                    by = c("Site" = "Site"))
rownames(bac_euk_veg$map_loaded) <- rn3
ncol(bac_euk_veg$map_loaded)
match("acacia_ataxacantha", names(bac_euk_veg$map_loaded))
veg <- bac_euk_veg$map_loaded[,match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
colnames(veg) == colnames(plantspcomp)
elli_veg <- filter_data(bac_euk_veg,
                        filter_cat = "subspecies",
                        keep_vals = "ellioti") # 28
schwein_veg <- filter_data(bac_euk_veg,
                           filter_cat = "subspecies",
                           keep_vals = "schweinfurthii") # 134
trog_veg <- filter_data(bac_euk_veg, 
                        filter_cat = "subspecies",
                        keep_vals = "troglodytes") # 35
verus_veg <- filter_data(bac_euk_veg, 
                         filter_cat = "subspecies",
                         keep_vals = "verus") # 312
veg.e <- elli_veg$map_loaded[,match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.s <- schwein_veg$map_loaded[,match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.t <- trog_veg$map_loaded[,match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.v <- verus_veg$map_loaded[,match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg_rel <- decostand(veg, method = "total")
veg_rel.e <- decostand(veg.e, method = "total")
veg_rel.s <- decostand(veg.s, method = "total")
veg_rel.t <- decostand(veg.t, method = "total")
veg_rel.v <- decostand(veg.v, method = "total")

bacteria.distance.v <- calc_dm(bac_euk_veg$data_loaded, method = "bray_sq_trans")
elli.bacteria.distance.v <- calc_dm(elli_veg$data_loaded, method = "bray_sq_trans")
schwein.bacteria.distance.v <- calc_dm(schwein_veg$data_loaded, method = "bray_sq_trans")
trog.bacteria.distance.v <- calc_dm(trog_veg$data_loaded, method = "bray_sq_trans")
verus.bacteria.distance.v <- calc_dm(verus_veg$data_loaded, method = "bray_sq_trans")

veg.distance <- vegdist(veg_rel, method = "bray")
elli.veg.distance <- vegdist(veg_rel.e, method = "bray")
schwein.veg.distance <- vegdist(veg_rel.s, method = "bray")
trog.veg.distance <- vegdist(veg_rel.t, method = "bray")
verus.veg.distance <- vegdist(veg_rel.v, method = "bray")

geog.v <- select(bac_euk_veg$map_loaded, long, lat)
geog.distance.v.mat <- rdist.earth(geog.v, miles = FALSE, R = 6371)
rownames(geog.distance.v.mat) <- colnames(geog.distance.v.mat) <- rownames(geog.distance.v.mat)
geog.distance.v <- as.dist(geog.distance.v.mat)
elli.geog.v <- select(elli_veg$map_loaded, long, lat)
elli.geog.distance.v.mat <- rdist.earth(elli.geog.v, miles = FALSE, R = 6371)
rownames(elli.geog.distance.v.mat) <- colnames(elli.geog.distance.v.mat) <- rownames(elli.geog.distance.v.mat)
elli.geog.distance.v <- as.dist(elli.geog.distance.v.mat)
schwein.geog.v <- select(schwein_veg$map_loaded, long, lat)
schwein.geog.distance.v.mat <- rdist.earth(schwein.geog.v, miles = FALSE, R = 6371)
rownames(schwein.geog.distance.v.mat) <- colnames(schwein.geog.distance.v.mat) <- rownames(schwein.geog.distance.v.mat)
schwein.geog.distance.v <- as.dist(schwein.geog.distance.v.mat)
trog.geog.v <- select(trog_veg$map_loaded, long, lat)
trog.geog.distance.v.mat <- rdist.earth(trog.geog.v, miles = FALSE, R = 6371)
rownames(trog.geog.distance.v.mat) <- colnames(trog.geog.distance.v.mat) <- rownames(trog.geog.distance.v.mat)
trog.geog.distance.v <- as.dist(trog.geog.distance.v.mat)
verus.geog.v <- select(verus_veg$map_loaded, long, lat)
verus.geog.distance.v.mat <- rdist.earth(verus.geog.v, miles = FALSE, R = 6371)
rownames(verus.geog.distance.v.mat) <- colnames(verus.geog.distance.v.mat) <- rownames(verus.geog.distance.v.mat)
verus.geog.distance.v <- as.dist(verus.geog.distance.v.mat)

# Make Parasite distance matrices with veg data
parasites.v <- select(bac_euk_veg$map_loaded,
                      Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                      Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                      Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                      Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                      Trichomonadidae_8, Troglodytella_abrassarti)
parasites.v.distance <- vegdist(parasites.v, method = "jaccard")
parasites.v.distance <- dist.zeroes(parasites.v, parasites.v.distance)
elli.parasites.v <- select(elli_veg$map_loaded, 
                           Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                           Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                           Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                           Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                           Trichomonadidae_8, Troglodytella_abrassarti)
elli.parasites.v.distance <- vegdist(elli.parasites.v, method = "jaccard")
elli.parasites.v.distance <- dist.zeroes(elli.parasites.v, elli.parasites.v.distance)
schwein.parasites.v <- select(schwein_veg$map_loaded,
                              Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                              Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                              Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                              Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                              Trichomonadidae_8, Troglodytella_abrassarti)
schwein.parasites.v.distance <- vegdist(schwein.parasites.v, method = "jaccard")
trog.parasites.v <- select(trog_veg$map_loaded, 
                           Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                           Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                           Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                           Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                           Trichomonadidae_8, Troglodytella_abrassarti)
trog.parasites.v.distance <- vegdist(trog.parasites.v, method = "jaccard")
verus.parasites.v <- select(verus_veg$map_loaded, 
                            Blastocystis_59, Blepharocorys_uncinata, Chilomastix_mesnili,
                            Dientamoeba_fragilis, Entamoeba_hartmanni, Entamoeba_muris, 
                            Haemonchus_contortus,Strongyloides_fuelleborni, Tetratrichomonas_2, 
                            Tetratrichomonas_31, Trichomonadidae_12_51, Trichomonadidae_15, 
                            Trichomonadidae_8, Troglodytella_abrassarti)
verus.parasites.v.distance <- vegdist(verus.parasites.v, method = "jaccard")

# Sampling
range(bac_euk$map_loaded$collectionyear)
table(bac_euk$map_loaded$collectionmonth)
range(bac_euk$map_loaded$collectionday)



########################### Supplemental Analyses and Figures #################################
############################## _1S. Alpha Rarefaction #########################################
samplemeans <- matrix(NA, nrow = 80, ncol = ncol(ad))
samplemeans[,2] <- aggregate(ad[,4], list(ad[,2]), mean)[,1]
for (i in 4:ncol(samplemeans)) {
  samplemeans[,i] <- aggregate(ad[,i], list(ad[,2]), mean)[,2]
}
samplemeans <- as.data.frame(samplemeans)
colnames(samplemeans) <- colnames(ad)
samplemeans <- select(samplemeans, -X, -iteration)
samplemeans$sequences.per.sample <- as.factor(samplemeans$sequences.per.sample)
longdata <- melt(samplemeans)
longdata$sequences.per.sample <- as.character(longdata$sequences.per.sample) %>% as.numeric(longdata$sequences.per.sample)
samples_subspecies <- bac_euk$map_loaded %>% select(sampleID, subspecies)
longdata <- inner_join(longdata, samples_subspecies, by = c("variable" = "sampleID"))
longdata <- na.omit(longdata)
longdata_mse <- ddply(longdata, c("subspecies","sequences.per.sample"), summarise,
                      seqs_mean = mean(value, na.rm = T),
                      seqs_se = se(value))
pdf(file = "SupplementaryFigureS1.pdf", width = 6, height = 4)
ggplot(longdata_mse, aes(sequences.per.sample, seqs_mean)) +
  geom_errorbar(aes(ymin = seqs_mean-seqs_se, ymax = seqs_mean+seqs_se, colour=subspecies),
                width=50, size=.5) +
  geom_line(aes(colour = subspecies, group = subspecies), size = 1) +
  labs(x = "Sequences per Sample",
       y = "# Observed ASVs",
       colour = "Region") +
  scale_color_discrete(breaks = c("ellioti", "troglodytes", "schweinfurthii", "verus"),
                       labels = c("Nigeria-Cameroon", "Central", "East", "West")) +
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = c(0.85,0.175),
        legend.title = element_text(size = 10),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.spacing.y = unit(0.1, "cm"))
dev.off()

# Mean richness by region
aggregate(longdata_mse$seqs_mean, 
          by = list(longdata_mse$subspecies), max, na.rm = TRUE)



########################## _2S. Within individual versus among individuals ####################
# Within sites (and subspecies)
# Just for individuals with repeat sampling
repeat_bray <- calc_dm(bac_input_filt_rar_repeat$data_loaded, method = "bray_sq_trans")
repeat_bray_mat <- as.matrix(repeat_bray)
repeat_bray_mat_df <- as.data.frame(repeat_bray_mat)
repeat_bray_mat_df$sample <- rownames(repeat_bray_mat_df)
repeat_bray_mat_df_melt <- melt(repeat_bray_mat_df, id.vars = "sample")
# 178084, same as matrix (422*422)
repeat_bray_mat_df_melt$sample <- as.factor(repeat_bray_mat_df_melt$sample)
# Get rid of zeros (same sample pairs)
repeat_bray_mat_df_melt <- subset(repeat_bray_mat_df_melt, value != 0)
droplevels(repeat_bray_mat_df_melt)
# 177662. 178084-177662 = 422. Good (should be 1 zero for each sample)
# Get rid of duplications (because upper matrix triangle repeats lower matrix triangle)
repeat_bray_mat_df_melt <- repeat_bray_mat_df_melt[ !duplicated(repeat_bray_mat_df_melt$value), ]
# 88831. 177662/2 = 88831. Good (should have cut dataframe in half)
droplevels(repeat_bray_mat_df_melt)
# Now we need to add the individual data
# This must be done twice, once to match col1, once to match col2
samp_ind <- bac_input_filt_rar_repeat$map_loaded %>% select(sampleID, consensusID)
repeat_bray_mat_df_melt_join <- inner_join(repeat_bray_mat_df_melt, samp_ind, 
                                           by = c("sample" = "sampleID"))
# Now col2
repeat_bray_mat_df_melt_join <- inner_join(repeat_bray_mat_df_melt_join, samp_ind, 
                                           by = c("variable" = "sampleID"))
# Make new column indicating if comparison is within or between individuals
for (i in 1:nrow(repeat_bray_mat_df_melt_join)) {
  ifelse(repeat_bray_mat_df_melt_join$consensusID.x[i] == repeat_bray_mat_df_melt_join$consensusID.y[i],
         repeat_bray_mat_df_melt_join$comparison[i] <- "within",
         repeat_bray_mat_df_melt_join$comparison[i] <- "between")
}
repeat_bray_mat_df_melt_join$comparison <- as.factor(repeat_bray_mat_df_melt_join$comparison)
table(repeat_bray_mat_df_melt_join$comparison)
# Add subspecies and site
samp_sp_site <- bac_input_filt_rar_repeat$map_loaded %>% select(sampleID, subspecies, Site)
repeat_bray_mat_df_melt_join <- inner_join(repeat_bray_mat_df_melt_join, samp_sp_site, 
                                           by = c("sample" = "sampleID"))
repeat_bray_mat_df_melt_join <- inner_join(repeat_bray_mat_df_melt_join, samp_sp_site, 
                                           by = c("variable" = "sampleID"))
# Make new column for within and among only within the same subspecies and site
repeat_bray_mat_df_melt_join$samesubspecies <- NA
for (i in 1:nrow(repeat_bray_mat_df_melt_join)) {
  ifelse(repeat_bray_mat_df_melt_join$subspecies.x[i] == repeat_bray_mat_df_melt_join$subspecies.y[i],
         repeat_bray_mat_df_melt_join$samesubspecies[i] <- "yes",
         repeat_bray_mat_df_melt_join$samesubspecies[i] <- "no")
}
repeat_bray_mat_df_melt_join$samesite <- NA
for (i in 1:nrow(repeat_bray_mat_df_melt_join)) {
  ifelse(repeat_bray_mat_df_melt_join$Site.x[i] == repeat_bray_mat_df_melt_join$Site.y[i],
         repeat_bray_mat_df_melt_join$samesite[i] <- "yes",
         repeat_bray_mat_df_melt_join$samesite[i] <- "no")
}
repeat_bray_mat_df_melt_join <- droplevels(repeat_bray_mat_df_melt_join)
repeat_bc_final <- subset(repeat_bray_mat_df_melt_join,
                          samesubspecies == "yes" & samesite == "yes")
table(repeat_bc_final$comparison)

# Subspecies dfs
elli_bc <- subset(repeat_bc_final, subspecies.x == "ellioti")
schwein_bc <- subset(repeat_bc_final, subspecies.x == "schweinfurthii")
trog_bc <- subset(repeat_bc_final, subspecies.x == "troglodytes")
verus_bc <- subset(repeat_bc_final, subspecies.x == "verus")

# Stats (all significant, higher BC in among)
t.test(value ~ comparison, data = elli_bc)
t.test(value ~ comparison, data = schwein_bc)
t.test(value ~ comparison, data = trog_bc)
t.test(value ~ comparison, data = verus_bc)

# Graph
levels(repeat_bc_final$subspecies.x)
repeat_bc_final$subspecies.x <- factor(repeat_bc_final$subspecies.x,
                                       levels =c("verus","ellioti","troglodytes","schweinfurthii"))
label.df.2S <- data.frame(
  subspecies = c("ellioti","ellioti",
                 "schweinfurthii","schweinfurthii",
                 "troglodytes","troglodytes",
                 "verus","verus"),
  comparison = c("between","within","between","within","between","within","between","within"),
  Value = c(1,1,1,1,1,1,1,1),
  Sig = c("a","b","a","b","a","b","a","b"))
facet_names.2S <- c('ellioti' = "b) N-C",
                    'schweinfurthii' = "d) East",
                    'troglodytes' = "c) Central",
                    'verus' = "a) West")
pdf(file = "SupplementaryFigureS2.pdf", width = 7, height = 3)
ggplot(data = repeat_bc_final, aes(comparison, value)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.05) +
  geom_text(data = label.df.2S, aes(x = comparison, y = Value, label = Sig, group = NULL)) +
  labs(x = "Individual Comparison",
       y = "Bray-Curtis Dissimilarity") +
  facet_wrap(~ subspecies.x, ncol = 4, labeller = as_labeller(facet_names.2S)) +
  ylim(0,1) +
  theme(axis.title = element_text(face="bold", size = 16),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 14))
dev.off()



######################## _3S. Geography vs. Genetics, climate, diet ###########################
# Climate
mantel(sl.geography.distance, sl.climate.distance) # p = 0.001, r = 0.45
mantel(sl.elli.geography.distance, sl.elli.climate.distance) # NA
mantel(sl.schwein.geography.distance, sl.schwein.climate.distance) # p = 0.017, r = 0.61
mantel(sl.trog.geography.distance, sl.trog.climate.distance) # p = 0.008, r = 0.62
mantel(sl.verus.geography.distance, sl.verus.climate.distance) # p = 0.001, r = 0.70

sl.geography.distance.mat[upper.tri(sl.geography.distance.mat, diag = TRUE)] <- NA
sl.geography.distance.mat.long <- melt(sl.geography.distance.mat)
sl.geography.distance.mat.long <- na.omit(sl.geography.distance.mat.long)
sl.climate.distance.mat <- as.matrix(sl.climate.distance)
sl.climate.distance.mat[upper.tri(sl.climate.distance.mat, diag = TRUE)] <- NA
sl.climate.distance.mat.long <- melt(sl.climate.distance.mat)
sl.climate.distance.mat.long <- na.omit(sl.climate.distance.mat.long)
sl.all.mat.long <- sl.geography.distance.mat.long
names(sl.all.mat.long) <- c("Site1","Site2","geog.dist.km")
sl.all.mat.long$clim.dist.euc <- sl.climate.distance.mat.long$value

sl.schwein.geography.distance.mat[upper.tri(sl.schwein.geography.distance.mat, diag = TRUE)]<-NA
sl.schwein.geography.distance.mat.long <- melt(sl.schwein.geography.distance.mat)
sl.schwein.geography.distance.mat.long <- na.omit(sl.schwein.geography.distance.mat.long)
sl.schwein.climate.distance.mat <- as.matrix(sl.schwein.climate.distance)
sl.schwein.climate.distance.mat[upper.tri(sl.schwein.climate.distance.mat, diag = TRUE)] <- NA
sl.schwein.climate.distance.mat.long <- melt(sl.schwein.climate.distance.mat)
sl.schwein.climate.distance.mat.long <- na.omit(sl.schwein.climate.distance.mat.long)
sl.schwein.all.mat.long <- sl.schwein.geography.distance.mat.long
names(sl.schwein.all.mat.long) <- c("Site1","Site2","geog.dist.km")
sl.schwein.all.mat.long$clim.dist.euc <- sl.schwein.climate.distance.mat.long$value

sl.trog.geography.distance.mat[upper.tri(sl.trog.geography.distance.mat, diag = TRUE)] <- NA
sl.trog.geography.distance.mat.long <- melt(sl.trog.geography.distance.mat)
sl.trog.geography.distance.mat.long <- na.omit(sl.trog.geography.distance.mat.long)
sl.trog.climate.distance.mat <- as.matrix(sl.trog.climate.distance)
sl.trog.climate.distance.mat[upper.tri(sl.trog.climate.distance.mat, diag = TRUE)] <- NA
sl.trog.climate.distance.mat.long <- melt(sl.trog.climate.distance.mat)
sl.trog.climate.distance.mat.long <- na.omit(sl.trog.climate.distance.mat.long)
sl.trog.all.mat.long <- sl.trog.geography.distance.mat.long
names(sl.trog.all.mat.long) <- c("Site1","Site2","geog.dist.km")
sl.trog.all.mat.long$clim.dist.euc <- sl.trog.climate.distance.mat.long$value

sl.verus.geography.distance.mat[upper.tri(sl.verus.geography.distance.mat, diag = TRUE)] <- NA
sl.verus.geography.distance.mat.long <- melt(sl.verus.geography.distance.mat)
sl.verus.geography.distance.mat.long <- na.omit(sl.verus.geography.distance.mat.long)
sl.verus.climate.distance.mat <- as.matrix(sl.verus.climate.distance)
sl.verus.climate.distance.mat[upper.tri(sl.verus.climate.distance.mat, diag = TRUE)] <- NA
sl.verus.climate.distance.mat.long <- melt(sl.verus.climate.distance.mat)
sl.verus.climate.distance.mat.long <- na.omit(sl.verus.climate.distance.mat.long)
sl.verus.all.mat.long <- sl.verus.geography.distance.mat.long
names(sl.verus.all.mat.long) <- c("Site1","Site2","geog.dist.km")
sl.verus.all.mat.long$clim.dist.euc <- sl.verus.climate.distance.mat.long$value

# Diet
mantel(sl.geography.distance, sl.diet.distance) # p = 0.01, r = 0.17
mantel(sl.elli.geography.distance, sl.elli.diet.distance) # NA
mantel(sl.schwein.geography.distance, sl.schwein.diet.distance) # p = 0.08, r = 0.71
mantel(sl.trog.geography.distance, sl.trog.diet.distance) # p = 0.63, r = -0.29
mantel(sl.verus.geography.distance, sl.verus.diet.distance) # p = 0.001, r = 0.51

sl.diet.distance.mat <- as.matrix(sl.diet.distance)
sl.diet.distance.mat[upper.tri(sl.diet.distance.mat, diag = TRUE)] <- NA
sl.diet.distance.mat.long <- melt(sl.diet.distance.mat)
sl.diet.distance.mat.long <- na.omit(sl.diet.distance.mat.long)
sl.all.mat.long$diet.dist.jac <- sl.diet.distance.mat.long$value

sl.schwein.diet.distance.mat <- as.matrix(sl.schwein.diet.distance)
sl.schwein.diet.distance.mat[upper.tri(sl.schwein.diet.distance.mat, diag = TRUE)] <- NA
sl.schwein.diet.distance.mat.long <- melt(sl.schwein.diet.distance.mat)
sl.schwein.diet.distance.mat.long <- na.omit(sl.schwein.diet.distance.mat.long)
sl.schwein.all.mat.long$diet.dist.jac <- sl.schwein.diet.distance.mat.long$value

sl.trog.diet.distance.mat <- as.matrix(sl.trog.diet.distance)
sl.trog.diet.distance.mat[upper.tri(sl.trog.diet.distance.mat, diag = TRUE)] <- NA
sl.trog.diet.distance.mat.long <- melt(sl.trog.diet.distance.mat)
sl.trog.diet.distance.mat.long <- na.omit(sl.trog.diet.distance.mat.long)
sl.trog.all.mat.long$diet.dist.jac <- sl.trog.diet.distance.mat.long$value

sl.verus.diet.distance.mat <- as.matrix(sl.verus.diet.distance)
sl.verus.diet.distance.mat[upper.tri(sl.verus.diet.distance.mat, diag = TRUE)] <- NA
sl.verus.diet.distance.mat.long <- melt(sl.verus.diet.distance.mat)
sl.verus.diet.distance.mat.long <- na.omit(sl.verus.diet.distance.mat.long)
sl.verus.all.mat.long$diet.dist.jac <- sl.verus.diet.distance.mat.long$value

# Genetics
gen_mat_df$site <- as.factor(rownames(gen_mat_df))
gen_mat_df_melt <- melt(gen_mat_df, id.vars = "site")
levels(gen_mat_df_melt$site)
levels(gen_mat_df_melt$variable)
gen_mat_df_melt$site <- str_replace_all(gen_mat_df_melt$site, '-', '.')
gen_mat_df_melt$site <- str_replace_all(gen_mat_df_melt$site, '_', '.')
gen_mat_df_melt$variable <- str_replace_all(gen_mat_df_melt$variable, '_', '.')
gen_mat_df_melt$site <- as.factor(gen_mat_df_melt$site)
gen_mat_df_melt$variable <- as.factor(gen_mat_df_melt$variable)
levels(gen_mat_df_melt$site)
levels(gen_mat_df_melt$variable)
gen_mat_df_melt <- gen_mat_df_melt[!is.na(gen_mat_df_melt$value),]
long_lat_site_sp <- select(bac_euk$map_loaded,
                        long, lat, Site, subspecies)
long_lat_site_sp <- long_lat_site_sp[!duplicated(long_lat_site_sp$lat),]
levels(long_lat_site_sp$Site)
long_lat_site_sp$Site <- str_replace_all(long_lat_site_sp$Site, '-', '.')
long_lat_site_sp$Site <- str_replace_all(long_lat_site_sp$Site, '_', '.')
long_lat_site_sp$Site <- as.factor(long_lat_site_sp$Site)
levels(long_lat_site_sp$Site)
geo_gen <- inner_join(gen_mat_df_melt, long_lat_site_sp, by = c("site" = "Site"))
geo_gen <- inner_join(geo_gen, long_lat_site_sp, by = c("variable" = "Site"))
# 276 = (24*23)/2 good! (24 matching sites in the genetic dataset)
# Distance in degrees
for (i in 1:nrow(geo_gen)) {
  y <- as.vector(c(geo_gen$long.x[i], geo_gen$lat.x[i]))
  z <- as.vector(c(geo_gen$long.y[i], geo_gen$lat.y[i]))
  geo_gen$geog.dist.deg[i] <- EuclideanDistance(y,z)
}

# Distance in km
geo_gen$geog.dist.km = mapply(lat_a=geo_gen$lat.x,
                              lon_a=geo_gen$long.x, 
                              lat_b=geo_gen$lat.y, 
                              lon_b=geo_gen$long.y, 
                              FUN = dist_geo)
max(geo_gen$geog.dist.km)
5300/2
min(geo_gen$geog.dist.km)

# Correlations and labels for all, and each subsp.
cor.test(geo_gen$geog.dist.km, geo_gen$value, 
         method = "pearson") # r = 0.65, p < 0.0001
elli_geo_gen <- subset(geo_gen, 
                       subspecies.x == "ellioti" & subspecies.y == "ellioti") #0
schwein_geo_gen <- subset(geo_gen, 
                          subspecies.x == "schweinfurthii"&subspecies.y == "schweinfurthii") #10
cor.test(schwein_geo_gen$geog.dist.km, schwein_geo_gen$value, 
         method = "pearson") # r = 0.20, p = 0.58
trog_geo_gen <- subset(geo_gen, 
                       subspecies.x == "troglodytes" & subspecies.y == "troglodytes") #10
cor.test(trog_geo_gen$geog.dist.km, trog_geo_gen$value, 
         method = "pearson") # r = 0.66, p = 0.04
verus_geo_gen <- subset(geo_gen, 
                       subspecies.x == "verus" & subspecies.y == "verus") #78
cor.test(verus_geo_gen$geog.dist.km, verus_geo_gen$value, 
         method = "pearson") # r = 0.19, p = 0.10

# Labels
label.df.3S.1 <- data.frame(label=c("n = 24, r = 0.64, p < 0.001"),
                            x = 3250, 
                            y = 0.01)
label.df.3S.2 <- data.frame(label=c("n = 5, r = 0.20, p = 0.58"),
                            x = 3250, 
                            y = 0.01)
label.df.3S.3 <- data.frame(label=c("n = 5, r = 0.66, p = 0.04"),
                            x = 3250, 
                            y = 0.01)
label.df.3S.4 <- data.frame(label=c("n = 13, r = 0.19, p = 0.10"),
                            x = 3250, 
                            y = 0.01)
label.df.3S.5 <- data.frame(label=c("n = 29, r = 0.45, p = 0.001"),
                            x = 3250, 
                            y = 0.1)
label.df.3S.6 <- data.frame(label=c("n = 5, r = 0.61, p = 0.001"),
                            x = 3250, 
                            y = 0.1)
label.df.3S.7 <- data.frame(label=c("n = 5, r = 0.62, p = 0.008"),
                            x = 3250, 
                            y = 0.1)
label.df.3S.8 <- data.frame(label=c("n = 17, r = 0.70, p = 0.001"),
                            x = 3250, 
                            y = 0.1)
label.df.3S.9 <- data.frame(label=c("n = 29, r = 0.17, p = 0.01"),
                            x = 3250, 
                            y = 0.01)
label.df.3S.10 <- data.frame(label=c("n = 5, r = 0.71, p = 0.08"),
                             x = 3250, 
                             y = 0.01)
label.df.3S.11 <- data.frame(label=c("n = 5, r = -0.29, p = 0.63"),
                             x = 3250, 
                             y = 0.01)
label.df.3S.12 <- data.frame(label=c("n = 17, r = 0.51, p = 0.001"),
                             x = 3250, 
                             y = 0.01)

# Graphs. Don't Graph ellioti because meaningless because only 2 sites!
g1 <- ggplot(geo_gen, aes(geog.dist.km, value)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.1, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 0.78) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g2 <- ggplot(schwein_geo_gen, aes(geog.dist.km, value)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_text(data = label.df.3S.2, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 0.78) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g3 <- ggplot(trog_geo_gen, aes(geog.dist.km, value)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.3, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 0.78) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g4 <- ggplot(verus_geo_gen, aes(geog.dist.km, value)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_text(data = label.df.3S.4, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 0.78) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g5 <- ggplot(sl.all.mat.long, aes(geog.dist.km, clim.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.5, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 6) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g6 <- ggplot(sl.schwein.all.mat.long, aes(geog.dist.km, clim.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.6, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 6) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g7 <- ggplot(sl.trog.all.mat.long, aes(geog.dist.km, clim.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.7, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 6) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g8 <- ggplot(sl.verus.all.mat.long, aes(geog.dist.km, clim.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.8, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 6) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g9 <- ggplot(sl.all.mat.long, aes(geog.dist.km, diet.dist.jac)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.9, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g10 <- ggplot(sl.schwein.all.mat.long, aes(geog.dist.km, diet.dist.jac)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_text(data = label.df.3S.10, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g11 <- ggplot(sl.trog.all.mat.long, aes(geog.dist.km, diet.dist.jac)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_text(data = label.df.3S.11, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g12 <- ggplot(sl.verus.all.mat.long, aes(geog.dist.km, diet.dist.jac)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.12, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

# Final figure will be made in powerpoint.
pdf(file = "forPPTS3.pdf", width = 10, height = 7)
plot_grid(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12, nrow = 3, align = "hv")
dev.off()

#### __Update: Add vegetation and remove diet #####
# Need to make site level geog and veg matrices
# Geog
geog.v.sl <- geog.v[ !duplicated(geog.v$long), ]
geog.distance.v.sl.mat <- rdist.earth(geog.v.sl, miles = FALSE, R = 6371)
rownames(geog.distance.v.sl.mat) <- colnames(geog.distance.v.sl.mat) <- rownames(geog.distance.v.sl.mat)
geog.distance.v.sl <- as.dist(geog.distance.v.sl.mat)
elli.geog.v.sl <- elli.geog.v[ !duplicated(elli.geog.v$long), ]
elli.geog.distance.v.sl.mat <- rdist.earth(elli.geog.v.sl, miles = FALSE, R = 6371)
rownames(elli.geog.distance.v.sl.mat) <- colnames(elli.geog.distance.v.sl.mat) <- rownames(elli.geog.distance.v.sl.mat)
elli.geog.distance.v.sl <- as.dist(elli.geog.distance.v.sl.mat)
schwein.geog.v.sl <- schwein.geog.v[ !duplicated(schwein.geog.v$long), ]
schwein.geog.distance.v.sl.mat <- rdist.earth(schwein.geog.v.sl, miles = FALSE, R = 6371)
rownames(schwein.geog.distance.v.sl.mat) <- colnames(schwein.geog.distance.v.sl.mat) <- rownames(schwein.geog.distance.v.sl.mat)
schwein.geog.distance.v.sl <- as.dist(schwein.geog.distance.v.sl.mat)
trog.geog.v.sl <- trog.geog.v[ !duplicated(trog.geog.v$long), ]
trog.geog.distance.v.sl.mat <- rdist.earth(trog.geog.v.sl, miles = FALSE, R = 6371)
rownames(trog.geog.distance.v.sl.mat) <- colnames(trog.geog.distance.v.sl.mat) <- rownames(trog.geog.distance.v.sl.mat)
trog.geog.distance.v.sl <- as.dist(trog.geog.distance.v.sl.mat)
verus.geog.v.sl <- verus.geog.v[ !duplicated(verus.geog.v$long), ]
verus.geog.distance.v.sl.mat <- rdist.earth(verus.geog.v.sl, miles = FALSE, R = 6371)
rownames(verus.geog.distance.v.sl.mat) <- colnames(verus.geog.distance.v.sl.mat) <- rownames(verus.geog.distance.v.sl.mat)
verus.geog.distance.v.sl <- as.dist(verus.geog.distance.v.sl.mat)
# Veg
veg.sl <- bac_euk_veg$map_loaded[!duplicated(bac_euk_veg$map_loaded$Site),match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.e.sl <- elli_veg$map_loaded[!duplicated(elli_veg$map_loaded$Site),match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.s.sl <- schwein_veg$map_loaded[!duplicated(schwein_veg$map_loaded$Site),match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.t.sl <- trog_veg$map_loaded[!duplicated(trog_veg$map_loaded$Site),match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg.v.sl <- verus_veg$map_loaded[!duplicated(verus_veg$map_loaded$Site),match("acacia_ataxacantha", names(bac_euk_veg$map_loaded)):ncol(bac_euk_veg$map_loaded)]
veg_rel.sl <- decostand(veg.sl, method = "total")
veg_rel.e.sl <- decostand(veg.e.sl, method = "total")
veg_rel.s.sl <- decostand(veg.s.sl, method = "total")
veg_rel.t.sl <- decostand(veg.t.sl, method = "total")
veg_rel.v.sl <- decostand(veg.v.sl, method = "total")
veg.distance.sl <- vegdist(veg_rel.sl, method = "bray")
elli.veg.distance.sl <- vegdist(veg_rel.e.sl, method = "bray")
schwein.veg.distance.sl <- vegdist(veg_rel.s.sl, method = "bray")
trog.veg.distance.sl <- vegdist(veg_rel.t.sl, method = "bray")
verus.veg.distance.sl <- vegdist(veg_rel.v.sl, method = "bray")

mantel(geog.distance.v.sl, veg.distance.sl) # p = 0.001, r = 0.52
mantel(elli.geog.distance.v.sl, elli.veg.distance.sl) # NA
mantel(schwein.geog.distance.v.sl, schwein.veg.distance.sl) # p = 1, r = -0.99
mantel(trog.geog.distance.v.sl, trog.veg.distance.sl) # p = 0.001, r = 0.82
mantel(verus.geog.distance.v.sl, verus.veg.distance.sl) # p = 0.001, r = 0.78

geog.distance.v.mat <- as.matrix(geog.distance.v.sl)
geog.distance.v.mat[upper.tri(geog.distance.v.mat, diag = TRUE)] <- NA
geog.distance.v.mat.long <- melt(geog.distance.v.mat)
geog.distance.v.mat.long <- na.omit(geog.distance.v.mat.long)
veg.distance.mat <- as.matrix(veg.distance.sl)
veg.distance.mat[upper.tri(veg.distance.mat, diag = TRUE)] <- NA
veg.distance.mat.long <- melt(veg.distance.mat)
veg.distance.mat.long <- na.omit(veg.distance.mat.long)
all.mat.long.v <- geog.distance.v.mat.long
names(all.mat.long.v) <- c("Site1","Site2","geog.dist.km")
all.mat.long.v$veg.dist.euc <- veg.distance.mat.long$value

schwein.geog.distance.v.mat <- as.matrix(schwein.geog.distance.v.sl)
schwein.geog.distance.v.mat[upper.tri(schwein.geog.distance.v.mat, diag = TRUE)]<-NA
schwein.geog.distance.v.mat.long <- melt(schwein.geog.distance.v.mat)
schwein.geog.distance.v.mat.long <- na.omit(schwein.geog.distance.v.mat.long)
schwein.veg.distance.mat <- as.matrix(schwein.veg.distance.sl)
schwein.veg.distance.mat[upper.tri(schwein.veg.distance.mat, diag = TRUE)] <- NA
schwein.veg.distance.mat.long <- melt(schwein.veg.distance.mat)
schwein.veg.distance.mat.long <- na.omit(schwein.veg.distance.mat.long)
schwein.all.mat.long.v <- schwein.geog.distance.v.mat.long
names(schwein.all.mat.long.v) <- c("Site1","Site2","geog.dist.km")
schwein.all.mat.long.v$veg.dist.euc <- schwein.veg.distance.mat.long$value

trog.geog.distance.v.mat <- as.matrix(trog.geog.distance.v.sl)
trog.geog.distance.v.mat[upper.tri(trog.geog.distance.v.mat, diag = TRUE)] <- NA
trog.geog.distance.v.mat.long <- melt(trog.geog.distance.v.mat)
trog.geog.distance.v.mat.long <- na.omit(trog.geog.distance.v.mat.long)
trog.veg.distance.mat <- as.matrix(trog.veg.distance.sl)
trog.veg.distance.mat[upper.tri(trog.veg.distance.mat, diag = TRUE)] <- NA
trog.veg.distance.mat.long <- melt(trog.veg.distance.mat)
trog.veg.distance.mat.long <- na.omit(trog.veg.distance.mat.long)
trog.all.mat.long.v <- trog.geog.distance.v.mat.long
names(trog.all.mat.long.v) <- c("Site1","Site2","geog.dist.km")
trog.all.mat.long.v$veg.dist.euc <- trog.veg.distance.mat.long$value

verus.geog.distance.v.mat <- as.matrix(verus.geog.distance.v.sl)
verus.geog.distance.v.mat[upper.tri(verus.geog.distance.v.mat, diag = TRUE)] <- NA
verus.geog.distance.v.mat.long <- melt(verus.geog.distance.v.mat)
verus.geog.distance.v.mat.long <- na.omit(verus.geog.distance.v.mat.long)
verus.veg.distance.mat <- as.matrix(verus.veg.distance.sl)
verus.veg.distance.mat[upper.tri(verus.veg.distance.mat, diag = TRUE)] <- NA
verus.veg.distance.mat.long <- melt(verus.veg.distance.mat)
verus.veg.distance.mat.long <- na.omit(verus.veg.distance.mat.long)
verus.all.mat.long.v <- verus.geog.distance.v.mat.long
names(verus.all.mat.long.v) <- c("Site1","Site2","geog.dist.km")
verus.all.mat.long.v$veg.dist.euc <- verus.veg.distance.mat.long$value

label.df.3S.13 <- data.frame(label=c("n = 27, r = 0.54, p = 0.001"),
                             x = 3250, 
                             y = 0.01)
label.df.3S.14 <- data.frame(label=c("n = 5, r = 0.81, p = 0.001"),
                             x = 3250, 
                             y = 0.01)
label.df.3S.15 <- data.frame(label=c("n = 3, r = -0.99, p = 1"),
                             x = 3250, 
                             y = 0.01)
label.df.3S.16 <- data.frame(label=c("n = 17, r = 0.78, p = 0.001"),
                             x = 3250, 
                             y = 0.01)

g13 <- ggplot(all.mat.long.v, aes(geog.dist.km, veg.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.13, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g14 <- ggplot(schwein.all.mat.long.v, aes(geog.dist.km, veg.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.14, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g15 <- ggplot(trog.all.mat.long.v, aes(geog.dist.km, veg.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_text(data = label.df.3S.15, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

g16 <- ggplot(verus.all.mat.long.v, aes(geog.dist.km, veg.dist.euc)) +
  geom_point(size = 3, alpha = 0.1) +
  geom_smooth() +
  geom_text(data = label.df.3S.16, aes(x = x, y = y, label = label), size = 4) +
  xlim(0, 5300) +
  ylim(0, 1) +
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank())

# Final figure will be made in powerpoint.
pdf(file = "forPPTS3v.pdf", width = 10, height = 5)
plot_grid(g1,g4,g3,g2,
          g5,g8,g7,g6,
          g13,g16,g15,g14, nrow = 3, align = "hv")
dev.off()



############################# _4S. Climate ####################################################
# Climate Distance and Bray-Curtis and Jaccard
# To combine with PCoA for Supplementary Climate Figure
# After running section 5S and 5, just take the climate parts and combine.
climpara <- subset(figS5comb, dataset == "Climate")
climprok <- subset(fig5comb, dataset == "Climate")
names(climpara)[3] <- "dissim"
names(climprok)[3] <- "dissim"
climpara$taxon <- "Parasites"
climprok$taxon <- "Prokaryotes"
climppt2 <- rbind(climprok,climpara)
climppt2$taxon <- factor(climppt2$taxon, levels = c("Prokaryotes","Parasites"))
(max(bacteria.distance.df.long$Climate) - min(bacteria.distance.df.long$Climate))*0.5
max(climppt2$dist)/2
label.df.clim <- data.frame(subspecies = c("all","all",
                                           "schweinfurthii","schweinfurthii",
                                           "troglodytes","troglodytes",
                                           "verus","verus"),
                            taxon = c("Prokaryotes","Parasites",
                                      "Prokaryotes","Parasites",
                                      "Prokaryotes","Parasites",
                                      "Prokaryotes","Parasites"),
                            x = c(3.16,3.16,3.16,3.16,
                                  3.16,3.16,3.16,3.16),
                            y = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08),
                            label = c("n = 560, r = 0.40, p = 0.001",
                                      "n = 560, r = 0.26, p = 0.001",
                                      "n = 134, r = 0.58, p = 0.001",
                                      "n = 134, r = 0.07, p = 0.004",
                                      "n = 86, r = 0.53, p = 0.001",
                                      "n = 86, r = 0.28, p = 0.001",
                                      "n = 312, r = 0.24, p = 0.001",
                                      "n = 312, r = 0.11, p = 0.001"))
facet_names_S4.v <- c("Prokaryotes" = "Prokaryotes", "Parasites" = "Parasites",
                      "all" = "all", "schweinfurthii" = "East", 
                      "troglodytes" = "Central", "verus" = "West")
levels(climppt2$subspecies)
climppt2$subspecies <- factor(climppt2$subspecies,
                              levels = c("all", "verus", "troglodytes", "schweinfurthii"))
pdf(file = "forPPTclim2.pdf", width = 3.5, height = 5)
ggplot(data = climppt2, aes(dist, dissim)) +
  geom_point(size = 0.5, alpha = 0.01) +
  geom_smooth(method = lm, se = T, size = 0.5) +
  geom_text(data = label.df.clim, aes(x = x, y = y, label = label), size = 2.5) +
  labs(x = "Euclidean Climate Distance",
       y = "Dissimilarity") +
  facet_grid(subspecies ~ taxon, labeller = as_labeller(facet_names_S4.v)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))
dev.off()



########################### _5S. Jaccard: Geography, Vegetation ##############################
# Run section 5 first.
# 8 panels (2 variables * all, schwein, trog, verus)
# No elli because only 2 sites (Geo, Climate, Diet the same), same as within vs among sites
# Not done at site level, done on whole dataset!
# Mantel Tests
mantel(parasites.distance, geography.distance) # p = 0.001, r = 0.45
mantel(elli.parasites.distance, elli.geography.distance) # p = 0.001, r = 0.54 (same for all)
mantel(schwein.parasites.distance, schwein.geography.distance) # p = 0.004, r = 0.12
mantel(trog.parasites.distance, trog.geography.distance) # p = 0.001, r = 0.35
mantel(verus.parasites.distance, verus.geography.distance) # p = 0.001, r = 0.11

mantel(parasites.distance, climate.distance) # p = 0.001, r = 0.26
mantel(elli.parasites.distance, elli.climate.distance) # 0.001, 0.54
mantel(schwein.parasites.distance, schwein.climate.distance) # p = 0.004, r = 0.07
mantel(trog.parasites.distance, trog.climate.distance) # p = 0.001, r = 0.28
mantel(verus.parasites.distance, verus.climate.distance) # p = 0.001, r = 0.11
mantel.partial(parasites.distance, climate.distance, geography.distance) # p = 0.22, r = 0.01
mantel.partial(elli.parasites.distance, elli.climate.distance, elli.geography.distance) # NA
mantel.partial(schwein.parasites.distance, 
               schwein.climate.distance, 
               schwein.geography.distance) # p = 0.96, r = -0.04
mantel.partial(trog.parasites.distance, 
               trog.climate.distance, 
               trog.geography.distance) # p = 0.56, r = -0.01
mantel.partial(verus.parasites.distance, 
               verus.climate.distance, 
               verus.geography.distance) # p = 0.13, r = 0.03

mantel(parasites.distance, diet.distance) # p = 0.001, r = 0.18
mantel(elli.parasites.distance, elli.diet.distance) # p = 0.001, r = 0.54
mantel(schwein.parasites.distance, schwein.diet.distance) # p = 0.01, r = 0.09
mantel(trog.parasites.distance, trog.diet.distance) # p = 0.009, r = 0.11
mantel(verus.parasites.distance, verus.diet.distance) # p = 0.001, r = 0.15
mantel.partial(parasites.distance, diet.distance, geography.distance) # p = 0.001, r = 0.06
mantel.partial(elli.parasites.distance, elli.diet.distance, elli.geography.distance) # NA
mantel.partial(schwein.parasites.distance, 
               schwein.diet.distance, 
               schwein.geography.distance) # p = 0.92, r = -0.02
mantel.partial(trog.parasites.distance, 
               trog.diet.distance, 
               trog.geography.distance) # p = 0.90, r = -0.07
mantel.partial(verus.parasites.distance, 
               verus.diet.distance, 
               verus.geography.distance) # p = 0.001, r = 0.10

# Dataframes (use same ones as in section 5. Must run section 5. first!)
# All
parasites.distance.mat <- as.matrix(parasites.distance)
parasites.distance.mat[upper.tri(parasites.distance.mat, diag = TRUE)] <- NA
parasites.distance.df <- as.data.frame(parasites.distance.mat)
parasites.distance.df$sampleID <- rownames(parasites.distance.df)
parasites.distance.df.long <- melt(parasites.distance.df, id.vars = "sampleID")
parasites.distance.df.long <- na.omit(parasites.distance.df.long)
names(parasites.distance.df.long) <- c("sampleID","variable","Jac")
parasites.distance.df.long$subspecies <- "all"
parasites.distance.df.long$Geography <- geography.distance.df.long$value
parasites.distance.df.long$Climate <- climate.distance.df.long$value
parasites.distance.df.long$Diet <- diet.distance.df.long$value
figS5all <- melt(parasites.distance.df.long,
                id.vars = c("sampleID","variable","Jac","subspecies"),
                measure.vars = c("Geography","Climate","Diet"))
names(figS5all) <- c("sampleID","variable","Jac","subspecies","dataset","dist")

# Schwein
schwein.parasites.distance.mat <- as.matrix(schwein.parasites.distance)
schwein.parasites.distance.mat[upper.tri(schwein.parasites.distance.mat, diag = TRUE)] <- NA
schwein.parasites.distance.df <- as.data.frame(schwein.parasites.distance.mat)
schwein.parasites.distance.df$sampleID <- rownames(schwein.parasites.distance.df)
schwein.parasites.distance.df.long <- melt(schwein.parasites.distance.df, id.vars = "sampleID")
schwein.parasites.distance.df.long <- na.omit(schwein.parasites.distance.df.long)
names(schwein.parasites.distance.df.long) <- c("sampleID","variable","Jac")
schwein.parasites.distance.df.long$subspecies <- "schweinfurthii"
schwein.parasites.distance.df.long$Geography <- schwein.geography.distance.df.long$value
schwein.parasites.distance.df.long$Climate <- schwein.climate.distance.df.long$value
schwein.parasites.distance.df.long$Diet <- schwein.diet.distance.df.long$value
figS5schwein <- melt(schwein.parasites.distance.df.long,
                    id.vars = c("sampleID","variable","Jac","subspecies"),
                    measure.vars = c("Geography","Climate","Diet"))
names(figS5schwein) <- c("sampleID","variable","Jac","subspecies","dataset","dist")

# Trog
trog.parasites.distance.mat <- as.matrix(trog.parasites.distance)
trog.parasites.distance.mat[upper.tri(trog.parasites.distance.mat, diag = TRUE)] <- NA
trog.parasites.distance.df <- as.data.frame(trog.parasites.distance.mat)
trog.parasites.distance.df$sampleID <- rownames(trog.parasites.distance.df)
trog.parasites.distance.df.long <- melt(trog.parasites.distance.df, id.vars = "sampleID")
trog.parasites.distance.df.long <- na.omit(trog.parasites.distance.df.long)
names(trog.parasites.distance.df.long) <- c("sampleID","variable","Jac")
trog.parasites.distance.df.long$subspecies <- "troglodytes"
trog.parasites.distance.df.long$Geography <- trog.geography.distance.df.long$value
trog.parasites.distance.df.long$Climate <- trog.climate.distance.df.long$value
trog.parasites.distance.df.long$Diet <- trog.diet.distance.df.long$value
figS5trog <- melt(trog.parasites.distance.df.long,
                 id.vars = c("sampleID","variable","Jac","subspecies"),
                 measure.vars = c("Geography","Climate","Diet"))
names(figS5trog) <- c("sampleID","variable","Jac","subspecies","dataset","dist")

# Verus
verus.parasites.distance.mat <- as.matrix(verus.parasites.distance)
verus.parasites.distance.mat[upper.tri(verus.parasites.distance.mat, diag = TRUE)] <- NA
verus.parasites.distance.df <- as.data.frame(verus.parasites.distance.mat)
verus.parasites.distance.df$sampleID <- rownames(verus.parasites.distance.df)
verus.parasites.distance.df.long <- melt(verus.parasites.distance.df, id.vars = "sampleID")
verus.parasites.distance.df.long <- na.omit(verus.parasites.distance.df.long)
names(verus.parasites.distance.df.long) <- c("sampleID","variable","Jac")
verus.parasites.distance.df.long$subspecies <- "verus"
verus.parasites.distance.df.long$Geography <- verus.geography.distance.df.long$value
verus.parasites.distance.df.long$Climate <- verus.climate.distance.df.long$value
verus.parasites.distance.df.long$Diet <- verus.diet.distance.df.long$value
figS5verus <- melt(verus.parasites.distance.df.long,
                  id.vars = c("sampleID","variable","Jac","subspecies"),
                  measure.vars = c("Geography","Climate","Diet"))
names(figS5verus) <- c("sampleID","variable","Jac","subspecies","dataset","dist")
figS5comb <- rbind(figS5all, figS5schwein, figS5trog, figS5verus)
figS5comb$subspecies <- as.factor(figS5comb$subspecies)

# Figure
(max(bacteria.distance.df.long$Geography) - min(bacteria.distance.df.long$Geography))*0.7
(max(bacteria.distance.df.long$Climate) - min(bacteria.distance.df.long$Climate))*0.75
(max(bacteria.distance.df.long$Diet) - min(bacteria.distance.df.long$Diet))*0.7
label.df.S5 <- data.frame(subspecies = c("all","all","all",
                                        "schweinfurthii","schweinfurthii","schweinfurthii",
                                        "troglodytes","troglodytes","troglodytes",
                                        "verus","verus","verus"),
                         dataset = c("Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet"),
                         x = c(3697.265, 4.40079, 0.7,
                               3697.265, 4.40079, 0.7,
                               3697.265, 4.40079, 0.7,
                               3697.265, 4.40079, 0.7),
                         y = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08),
                         label = c("n = 560, r = 0.45, p = 0.001",
                                   "n = 560, r = 0.26, p = 0.001",
                                   "n = 560, r = 0.18, p = 0.001",
                                   "n = 134, r = 0.12, p = 0.004",
                                   "n = 134, r = 0.07, p = 0.004",
                                   "n = 134, r = 0.09, p = 0.01",
                                   "n = 86, r = 0.35, p = 0.001",
                                   "n = 86, r = 0.28, p = 0.001",
                                   "n = 86, r = 0.11, p = 0.009",
                                   "n = 312, r = 0.11, p = 0.001",
                                   "n = 312, r = 0.11, p = 0.001",
                                   "n = 312, r = 0.15, p = 0.001"))
ggplot(data = figS5comb, aes(dist, Jac)) +
  geom_point(size = 0.5, alpha = 0.01) +
  geom_smooth(method = lm, se = T, size = 0.5) +
  geom_text(data = label.df.S5, aes(x = x, y = y, label = label), size = 2) +
  labs(x = "Distance",
       y = "Jaccard Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x") +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))



#### __Update: Add vegetation and remove climate and diet ####
mantel(veg.distance, parasites.v.distance) # r = 0.32, p = 0.001
mantel(elli.veg.distance, elli.parasites.v.distance) # r = 0.54, p = 0.001
mantel(schwein.veg.distance, schwein.parasites.v.distance) # r = 0.10, p = 0.001
mantel(trog.veg.distance, trog.parasites.v.distance) # r = 0.29, p = 0.001
mantel(verus.veg.distance, verus.parasites.v.distance) # r = 0.11, p = 0.001

mantel.partial(veg.distance, parasites.v.distance, geog.distance.v) # r = 0.08, p = 0.001
mantel.partial(elli.veg.distance,elli.parasites.v.distance,elli.geog.distance.v) # NA
mantel.partial(schwein.veg.distance,schwein.parasites.v.distance,schwein.geog.distance.v)#.002.45
mantel.partial(trog.veg.distance, trog.parasites.v.distance, trog.geog.distance.v) # 0.32, 0.001
mantel.partial(verus.veg.distance, verus.parasites.v.distance, verus.geog.distance.v) #0.04,0.033

# Assemble Figure Dataframe. Geography and Diet are the same. Use subset bac, veg instead of clim
parasites.v.distance.mat <- as.matrix(parasites.v.distance)
parasites.v.distance.mat[upper.tri(parasites.v.distance.mat, diag = TRUE)] <- NA
parasites.v.distance.df <- as.data.frame(parasites.v.distance.mat)
parasites.v.distance.df$sampleID <- rownames(parasites.v.distance.df)
parasites.v.distance.df.long <- melt(parasites.v.distance.df, id.vars = "sampleID")
parasites.v.distance.df.long <- na.omit(parasites.v.distance.df.long)
veg.distance.mat <- as.matrix(veg.distance)
veg.distance.mat[upper.tri(veg.distance.mat, diag = TRUE)] <- NA
veg.distance.df <- as.data.frame(veg.distance.mat)
veg.distance.df$sampleID <- rownames(veg.distance.df)
veg.distance.df.long <- melt(veg.distance.df, id.vars = "sampleID")
veg.distance.df.long <- na.omit(veg.distance.df.long)
names(parasites.v.distance.df.long) <- c("sampleID","variable","Jac")
parasites.v.distance.df.long$subspecies <- "all"
parasites.v.distance.df.long$Veg <- veg.distance.df.long$value
figS5all.vGD <- subset(figS5all, dataset != "Climate")
figS5all.vV <- melt(parasites.v.distance.df.long,
                    id.vars = c("sampleID","variable","Jac","subspecies"),
                    measure.vars = c("Veg"))
names(figS5all.vV) <- c("sampleID","variable","Jac","subspecies","dataset","dist")
figS5all.v <- rbind(figS5all.vGD, figS5all.vV)

# Schwein
schwein.parasites.v.distance.mat <- as.matrix(schwein.parasites.v.distance)
schwein.parasites.v.distance.mat[upper.tri(schwein.parasites.v.distance.mat, diag = TRUE)] <- NA
schwein.parasites.v.distance.df <- as.data.frame(schwein.parasites.v.distance.mat)
schwein.parasites.v.distance.df$sampleID <- rownames(schwein.parasites.v.distance.df)
schwein.parasites.v.distance.df.long<-melt(schwein.parasites.v.distance.df, id.vars = "sampleID")
schwein.parasites.v.distance.df.long <- na.omit(schwein.parasites.v.distance.df.long)
schwein.veg.distance.mat <- as.matrix(schwein.veg.distance)
schwein.veg.distance.mat[upper.tri(schwein.veg.distance.mat, diag = TRUE)] <- NA
schwein.veg.distance.df <- as.data.frame(schwein.veg.distance.mat)
schwein.veg.distance.df$sampleID <- rownames(schwein.veg.distance.df)
schwein.veg.distance.df.long <- melt(schwein.veg.distance.df, id.vars = "sampleID")
schwein.veg.distance.df.long <- na.omit(schwein.veg.distance.df.long)
names(schwein.parasites.v.distance.df.long) <- c("sampleID","variable","Jac")
schwein.parasites.v.distance.df.long$subspecies <- "schweinfurthii"
schwein.parasites.v.distance.df.long$Veg <- schwein.veg.distance.df.long$value
schwein.figS5all.vGD <- subset(figS5schwein, dataset != "Climate")
schwein.figS5all.vV <- melt(schwein.parasites.v.distance.df.long,
                            id.vars = c("sampleID","variable","Jac","subspecies"),
                            measure.vars = c("Veg"))
names(schwein.figS5all.vV) <- c("sampleID","variable","Jac","subspecies","dataset","dist")
schwein.figS5all.v <- rbind(schwein.figS5all.vGD, schwein.figS5all.vV)

# Trog
trog.parasites.v.distance.mat <- as.matrix(trog.parasites.v.distance)
trog.parasites.v.distance.mat[upper.tri(trog.parasites.v.distance.mat, diag = TRUE)] <- NA
trog.parasites.v.distance.df <- as.data.frame(trog.parasites.v.distance.mat)
trog.parasites.v.distance.df$sampleID <- rownames(trog.parasites.v.distance.df)
trog.parasites.v.distance.df.long <- melt(trog.parasites.v.distance.df, id.vars = "sampleID")
trog.parasites.v.distance.df.long <- na.omit(trog.parasites.v.distance.df.long)
trog.veg.distance.mat <- as.matrix(trog.veg.distance)
trog.veg.distance.mat[upper.tri(trog.veg.distance.mat, diag = TRUE)] <- NA
trog.veg.distance.df <- as.data.frame(trog.veg.distance.mat)
trog.veg.distance.df$sampleID <- rownames(trog.veg.distance.df)
trog.veg.distance.df.long <- melt(trog.veg.distance.df, id.vars = "sampleID")
trog.veg.distance.df.long <- na.omit(trog.veg.distance.df.long)
names(trog.parasites.v.distance.df.long) <- c("sampleID","variable","Jac")
trog.parasites.v.distance.df.long$subspecies <- "troglodytes"
trog.parasites.v.distance.df.long$Veg <- trog.veg.distance.df.long$value
trog.figS5all.vGD <- subset(figS5trog, dataset != "Climate")
trog.figS5all.vV <- melt(trog.parasites.v.distance.df.long,
                         id.vars = c("sampleID","variable","Jac","subspecies"),
                         measure.vars = c("Veg"))
names(trog.figS5all.vV) <- c("sampleID","variable","Jac","subspecies","dataset","dist")
trog.figS5all.v <- rbind(trog.figS5all.vGD, trog.figS5all.vV)

# Verus
verus.parasites.v.distance.mat <- as.matrix(verus.parasites.v.distance)
verus.parasites.v.distance.mat[upper.tri(verus.parasites.v.distance.mat, diag = TRUE)] <- NA
verus.parasites.v.distance.df <- as.data.frame(verus.parasites.v.distance.mat)
verus.parasites.v.distance.df$sampleID <- rownames(verus.parasites.v.distance.df)
verus.parasites.v.distance.df.long <- melt(verus.parasites.v.distance.df, id.vars = "sampleID")
verus.parasites.v.distance.df.long <- na.omit(verus.parasites.v.distance.df.long)
verus.veg.distance.mat <- as.matrix(verus.veg.distance)
verus.veg.distance.mat[upper.tri(verus.veg.distance.mat, diag = TRUE)] <- NA
verus.veg.distance.df <- as.data.frame(verus.veg.distance.mat)
verus.veg.distance.df$sampleID <- rownames(verus.veg.distance.df)
verus.veg.distance.df.long <- melt(verus.veg.distance.df, id.vars = "sampleID")
verus.veg.distance.df.long <- na.omit(verus.veg.distance.df.long)
names(verus.parasites.v.distance.df.long) <- c("sampleID","variable","Jac")
verus.parasites.v.distance.df.long$subspecies <- "verus"
verus.parasites.v.distance.df.long$Veg <- verus.veg.distance.df.long$value
verus.figS5all.vGD <- subset(figS5verus, dataset != "Climate")
verus.figS5all.vV <- melt(verus.parasites.v.distance.df.long,
                          id.vars = c("sampleID","variable","Jac","subspecies"),
                          measure.vars = c("Veg"))
names(verus.figS5all.vV) <- c("sampleID","variable","Jac","subspecies","dataset","dist")
verus.figS5all.v <- rbind(verus.figS5all.vGD, verus.figS5all.vV)
figS5comb.v <- rbind(figS5all.v, schwein.figS5all.v, trog.figS5all.v, verus.figS5all.v)
figS5comb.v$subspecies <- as.factor(figS5comb.v$subspecies)
levels(figS5comb.v$subspecies)
levels(figS5comb.v$dataset)
figS5comb.v$dataset <- factor(figS5comb.v$dataset, levels = c("Geography","Veg","Diet"))

# Figure
(max(parasites.v.distance.df.long$Veg) - min(parasites.v.distance.df.long$Veg))/2
# Update to remove diet
figS5comb.v.nd <- subset(figS5comb.v, dataset != "Diet")
label.df.S5.v <- data.frame(subspecies = c("all","all",
                                           "schweinfurthii","schweinfurthii",
                                           "troglodytes","troglodytes",
                                           "verus","verus"),
                            dataset = c("Geography","Veg",
                                        "Geography","Veg",
                                        "Geography","Veg",
                                        "Geography","Veg"),
                            x = c(3697.265, 0.5,
                                  3697.265, 0.5,
                                  3697.265, 0.5,
                                  3697.265, 0.5),
                            y = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08),
                            label = c("n = 560, r = 0.68, p = 0.001",
                                      "n = 509, r = 0.32, p = 0.001",
                                      "n = 134, r = 0.51, p = 0.001",
                                      "n = 134, r = 0.10, p = 0.001",
                                      "n = 86, r = 0.52, p = 0.001",
                                      "n = 35, r = 0.29, p = 0.001",
                                      "n = 312, r = 0.29, p = 0.001",
                                      "n = 312, r = 0.11, p = 0.001"))
facet_names_S5.v <- c("Geography" = "a) Geography", "Veg" = "b) Vegetation",
                      "all" = "all", "schweinfurthii" = "East", 
                      "troglodytes" = "Central", "verus" = "West")
levels(figS5comb.v.nd$subspecies)
figS5comb.v.nd$subspecies <- factor(figS5comb.v.nd$subspecies,
                                    levels = c("all", "verus", "troglodytes", "schweinfurthii"))
# pdf(file = "SupplementaryFigureS5.pdf", width = 5, height = 5)
png(file = "PNGs/SupplementaryFigureS5.png", width = 5, height = 5, units = "in", res = 300)
ggplot(data = figS5comb.v.nd, aes(dist, Jac)) +
  geom_point(size = 0.5, alpha = 0.01) +
  geom_smooth(data = subset(figS5comb.v.nd, dataset == "Geography"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(figS5comb.v.nd, dataset == "Veg"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_text(data = label.df.S5.v, aes(x = x, y = y, label = label), size = 2.5) +
  labs(x = "Distance",
       y = "Jaccard Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x", labeller = as_labeller(facet_names_S5.v)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10),
        panel.spacing.x = unit(1, "lines"))
dev.off()



################################ _6S. Region and Diet indicator taxa ########################
# Dataset (top 278, in 5% of samples and have 5 sequences)
input_rar_uniq_abund_ubiq <- filter_taxa_from_input(bac_euk,
                                                    taxa_IDs_to_keep = abund_ubiq278$asvID)
asv_abund_ubiq <- t(input_rar_uniq_abund_ubiq$data_loaded)
asv_abund_ubiq_df <- as.data.frame(asv_abund_ubiq)
diets_mp_taxonomy <- input_rar_uniq_abund_ubiq$taxonomy_loaded

#### __Region ####
### Indicators
set.seed(500)
all_mp <- multipatt(asv_abund_ubiq_df, 
                    input_rar_uniq_abund_ubiq$map_loaded$subspecies, 
                    func = "IndVal.g", control = how(nperm=999))
summary(all_mp)
mp_result_df <- all_mp$sign
mp_result_df <- na.omit(mp_result_df)
mp_result_df <- subset(mp_result_df, p.value < 0.05)
for (i in 1:nrow(mp_result_df)) {
  mp_result_df$numgroups[i] <- sum(mp_result_df$s.ellioti[i] +
                                     mp_result_df$s.schweinfurthii[i] +
                                     mp_result_df$s.troglodytes[i] +
                                     mp_result_df$s.verus[i])
}
mp_result_df <- subset(mp_result_df, numgroups == 1)
# These are the 82 ASVs assigned to one subspecies (7 elli, 15 schwein, 15 trog, 45 verus)
mp_result_df$asvID <- rownames(mp_result_df)
mp_result_df$indicator <- NA
mp_result_df$indicator <- paste0(mp_result_df$s.ellioti,
                                 mp_result_df$s.schweinfurthii,
                                 mp_result_df$s.troglodytes,
                                 mp_result_df$s.verus)
mp_result_df$indicator <- as.factor(mp_result_df$indicator)
mp_result_df$indicator <- revalue(mp_result_df$indicator, 
                                  c("1000" = "ellioti",
                                    "0100" = "schweinfurthii",
                                    "0010" = "troglodytes",
                                    "0001" = "verus"))
# Add taxonomy to the indicator results and export
taxonomy_df <- input_rar_uniq_abund_ubiq$taxonomy_loaded
taxonomy_df_tax <- taxonomy_df %>%
  unite(taxonomy1,taxonomy2,taxonomy3,taxonomy4,taxonomy5,taxonomy6,
        sep = ";")
mp_result_df <- inner_join(mp_result_df, taxonomy_df_tax, by = c("asvID" = "taxonomy7"))
# write.csv(mp_result_df, "subspecies_indicators_seeded.csv")


### Tree
# Make phyloseq object, Remove the 3 archaea, make circular tree
bac_euk278 <- filter_taxa_from_input(bac_euk, taxa_IDs_to_keep = tree$tip.label)
bac_euk278$taxonomy_loaded <- left_join(bac_euk278$taxonomy_loaded, mp_result_df,
                                        by = c("taxonomy7" = "asvID"))
row.names(bac_euk278$taxonomy_loaded) <- row.names(bac_euk278$data_loaded)
tax.mat <- as.matrix(bac_euk278$taxonomy_loaded)
otu_tab <- otu_table(bac_euk278$data_loaded, taxa_are_rows = TRUE)
tax_tab <- tax_table(tax.mat)
sam_tab <- sample_data(bac_euk278$map_loaded)
no_mix <- row.names(subset(bac_euk278$taxonomy_loaded, taxonomy7 != "ASV_256" & taxonomy7 != "ASV_253" & taxonomy7 != "ASV_184" & taxonomy7 != "ASV_209" & taxonomy7 != "ASV_53" & taxonomy7 != "ASV_174" & taxonomy7 != "ASV_115" & taxonomy7 != "ASV_187" & taxonomy7 != "ASV_24" & taxonomy7 != "ASV_27" & taxonomy7 != "ASV_60" & taxonomy7 != "ASV_1" & taxonomy7 != "ASV_23" & taxonomy7 != "ASV_78" & taxonomy7 != "ASV_157"))
ps <- phyloseq(otu_tab, tax_tab, sam_tab, tree)
ps <- prune_taxa(no_mix, ps)

# Presence and absence of being an indicator can be considered a trait and tested for phylo signal
# Need to do for each subspecies
library(picante)
# Make Dichotomous
phy <- multi2di(ps@phy_tree)
# Make dataframes and test for phylosignal
e.trait <- ps@tax_table %>%
  as.data.frame() %>%
  select(s.ellioti)
e.trait$s.ellioti <- as.numeric(as.character(e.trait$s.ellioti))
e.trait[is.na(e.trait)] <- 0
multiPhylosignal(e.trait, phy)
s.trait <- ps@tax_table %>%
  as.data.frame() %>%
  select(s.schweinfurthii)
s.trait$s.schweinfurthii <- as.numeric(as.character(s.trait$s.schweinfurthii))
s.trait[is.na(s.trait)] <- 0
multiPhylosignal(s.trait, phy)
t.trait <- ps@tax_table %>%
  as.data.frame() %>%
  select(s.troglodytes)
t.trait$s.troglodytes <- as.numeric(as.character(t.trait$s.troglodytes))
t.trait[is.na(t.trait)] <- 0
multiPhylosignal(t.trait, phy)
v.trait <- ps@tax_table %>%
  as.data.frame() %>%
  select(s.verus)
v.trait$s.verus <- as.numeric(as.character(v.trait$s.verus))
v.trait[is.na(v.trait)] <- 0
multiPhylosignal(v.trait, phy)
# None signficant. Can also try method specifically for categorical trait

# Get node numbers
ggtree(ps) + 
  geom_tippoint(aes(color=taxonomy2), size = 1.5) +
  geom_tiplab(aes(label = taxonomy7, hjust = -0.3), size = 1.5) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.2, size = 1.5) + 
  labs(color = "Phylum") +
  theme(legend.position = c(0.15, 0.65),
        legend.background = element_blank())

circ1 <- ggtree(ps, layout = "circular", branch.length = "none") + 
  geom_tippoint(aes(color=indicator), size = 1.5) +
  geom_cladelabel(node = 313, label = "Actinobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[1], hjust = "right", angle = 39) +
  geom_cladelabel(node = 333, label = "Bacteroidetes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[2], hjust = "left", angle = -58) +
  geom_cladelabel(node = 289, label = "Chloroflexi", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[3], hjust = "right", angle = -5) +
  geom_cladelabel(node = 285, label = "Cyanobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[4], hjust = "right", angle = -20) +
  geom_cladelabel(node = 423, label = "Fibrobacteres", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[5], hjust = "right", angle = 54) +
  geom_cladelabel(node = c(512, 440, 448), label = "Firmicutes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[6], hjust = "left", angle = 79) +
  geom_cladelabel(node = 291, label = "Kiritimatiellaeota", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[7], hjust = "right", angle = 4) +
  geom_cladelabel(node = 277, label = "Proteobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5,  color = hue_pal()(9)[8], hjust = "right", angle = -14) +
  geom_cladelabel(node = 304, label = "Spirochaetes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[9], hjust = "right", angle = 21) +
  scale_color_manual(values = c("#F8766D","#7CAE00","#00BFC4","#C77CFF","white"),
                     labels = c("N-C", "East", "Central", "West", "NA"),
                     na.translate = F) +
  labs(color = "Region") +
  geom_text(x = -45, y = -40, label = "a) Region Indicator Taxa", size = 6) +
  theme(legend.position = c(-0.05, 0.5),
        legend.background = element_blank(),
        plot.margin = unit(c(-0.5,-5,-0.5,-0.5), "cm")) +
  xlim(0,36)
circ1

#### __Diet ####
### Honey, Termites indicators within each region
elli_abund_ubiq <- filter_data(input_rar_uniq_abund_ubiq, 
                               filter_cat = "subspecies", 
                               keep_vals = "ellioti")
elli_abund_ubiq_df <- as.data.frame(t(elli_abund_ubiq$data_loaded))
schwein_abund_ubiq <- filter_data(input_rar_uniq_abund_ubiq, 
                                  filter_cat = "subspecies", 
                                  keep_vals = "schweinfurthii")
schwein_abund_ubiq_df <- as.data.frame(t(schwein_abund_ubiq$data_loaded))
trog_abund_ubiq <- filter_data(input_rar_uniq_abund_ubiq, 
                               filter_cat = "subspecies", 
                               keep_vals = "troglodytes")
trog_abund_ubiq_df <- as.data.frame(t(trog_abund_ubiq$data_loaded))
verus_abund_ubiq <- filter_data(input_rar_uniq_abund_ubiq, 
                                filter_cat = "subspecies", 
                                keep_vals = "verus")
verus_abund_ubiq_df <- as.data.frame(t(verus_abund_ubiq$data_loaded))

# Honey
set.seed(500)
elli_honey_mp <- multipatt(elli_abund_ubiq_df, 
                           elli_abund_ubiq$map_loaded$honey, func="IndVal.g", 
                           control = how(nperm=999))
summary(elli_honey_mp) # 12 no, 35 yes
elli_honey_results <- elli_honey_mp$sign
elli_honey_results$ASV <- rownames(elli_honey_results)
elli_honey_results <- subset(elli_honey_results, p.value <= 0.05)
set.seed(500)
schwein_honey_mp <- multipatt(schwein_abund_ubiq_df, 
                              schwein_abund_ubiq$map_loaded$honey, func="IndVal.g", 
                              control = how(nperm=999))
summary(schwein_honey_mp) # 35 no, 38 yes
schwein_honey_results <- schwein_honey_mp$sign
schwein_honey_results$ASV <- rownames(schwein_honey_results)
schwein_honey_results <- subset(schwein_honey_results, p.value <= 0.05)
set.seed(500)
trog_honey_mp <- multipatt(trog_abund_ubiq_df, 
                           trog_abund_ubiq$map_loaded$honey, func="IndVal.g", 
                           control = how(nperm=999))
summary(trog_honey_mp) # 28 no, 11 yes
trog_honey_results <- trog_honey_mp$sign
trog_honey_results$ASV <- rownames(trog_honey_results)
trog_honey_results <- subset(trog_honey_results, p.value <= 0.05)
set.seed(500)
verus_honey_mp <- multipatt(verus_abund_ubiq_df, 
                            verus_abund_ubiq$map_loaded$honey, func="IndVal.g", 
                            control = how(nperm=999))
summary(verus_honey_mp) # 3 no, 0 yes
verus_honey_results <- verus_honey_mp$sign
verus_honey_results$ASV <- rownames(verus_honey_results)
verus_honey_results <- subset(verus_honey_results, p.value <= 0.05)
mergeCols <- c("ASV","ASV")
honey_results <- elli_honey_results %>% 
  full_join(schwein_honey_results, by = mergeCols) %>%
  full_join(trog_honey_results, by = mergeCols) %>%
  full_join(verus_honey_results, by = mergeCols)
honey_results$numgroups_yes <- rowSums(honey_results[, c(2, 8, 13, 18)], na.rm = TRUE)
max(honey_results$numgroups_yes) # none with all 4
honey_results$numgroups_no <- rowSums(honey_results[, c(1, 7, 12, 17)], na.rm = TRUE)
max(honey_results$numgroups_no) # none with all 4
honey_venn <- honey_results %>% select(ASV, s.1.x, s.1.y, s.1.x.x, s.1.y.y)
colnames(honey_venn) <- c("ASV", "ellioti", "schweinfurthii", "troglodytes", "verus")
honey_venn[is.na(honey_venn)] <- 0
honey_venn[honey_venn == 0] <- "FALSE"
honey_venn[honey_venn == 1] <- "TRUE"
honey_venn_wTax <- inner_join(honey_venn, diets_mp_taxonomy, by = c("ASV" = "taxonomy7"))

# Termites
set.seed(500)
elli_termites_mp <- multipatt(elli_abund_ubiq_df, 
                           elli_abund_ubiq$map_loaded$termites, func="IndVal.g", 
                           control = how(nperm=999))
summary(elli_termites_mp) # 35 no, 12 yes
elli_termites_results <- elli_termites_mp$sign
elli_termites_results$ASV <- rownames(elli_termites_results)
elli_termites_results <- subset(elli_termites_results, p.value <= 0.05)
set.seed(500)
schwein_termites_mp <- multipatt(schwein_abund_ubiq_df, 
                              schwein_abund_ubiq$map_loaded$termites, func="IndVal.g", 
                              control = how(nperm=999))
summary(schwein_termites_mp) # 37 no, 21 yes
schwein_termites_results <- schwein_termites_mp$sign
schwein_termites_results$ASV <- rownames(schwein_termites_results)
schwein_termites_results <- subset(schwein_termites_results, p.value <= 0.05)
set.seed(500)
trog_termites_mp <- multipatt(trog_abund_ubiq_df, 
                           trog_abund_ubiq$map_loaded$termites, func="IndVal.g", 
                           control = how(nperm=999))
summary(trog_termites_mp) # 27 no, 26 yes
trog_termites_results <- trog_termites_mp$sign
trog_termites_results$ASV <- rownames(trog_termites_results)
trog_termites_results <- subset(trog_termites_results, p.value <= 0.05)
set.seed(500)
verus_termites_mp <- multipatt(verus_abund_ubiq_df, 
                            verus_abund_ubiq$map_loaded$termites, func="IndVal.g", 
                            control = how(nperm=999))
summary(verus_termites_mp) # 7 no, 2 yes
verus_termites_results <- verus_termites_mp$sign
verus_termites_results$ASV <- rownames(verus_termites_results)
verus_termites_results <- subset(verus_termites_results, p.value <= 0.05)
mergeCols <- c("ASV","ASV")
termites_results <- elli_termites_results %>% 
  full_join(schwein_termites_results, by = mergeCols) %>%
  full_join(trog_termites_results, by = mergeCols) %>%
  full_join(verus_termites_results, by = mergeCols)
termites_results$numgroups_yes <- rowSums(termites_results[, c(2, 8, 13, 18)], na.rm = TRUE)
max(termites_results$numgroups_yes) # none with all 4
termites_results$numgroups_no <- rowSums(termites_results[, c(1, 7, 12, 17)], na.rm = TRUE)
max(termites_results$numgroups_no) # none with all 4
termites_venn <- termites_results %>% select(ASV, s.1.x, s.1.y, s.1.x.x, s.1.y.y)
colnames(termites_venn) <- c("ASV", "ellioti", "schweinfurthii", "troglodytes", "verus")
termites_venn[is.na(termites_venn)] <- 0
termites_venn[termites_venn == 0] <- "FALSE"
termites_venn[termites_venn == 1] <- "TRUE"
termites_venn_wTax <- inner_join(termites_venn, diets_mp_taxonomy, by = c("ASV" = "taxonomy7"))

### Tree
abund_ubiq278_data <- left_join(abund_ubiq278, honey_venn_wTax[,1:5], by = c("asvID" = "ASV"))
abund_ubiq278_data <- left_join(abund_ubiq278_data, termites_venn_wTax[,1:5], 
                                by = c("asvID" = "ASV"))
abund_ubiq278_data[abund_ubiq278_data == TRUE] <- 1
abund_ubiq278_data[is.na(abund_ubiq278_data)] <- 0
abund_ubiq278_data[abund_ubiq278_data == FALSE] <- 0
abund_ubiq278_taxa <- left_join(abund_ubiq278, bac_euk$taxonomy_loaded, 
                                by = c("asvID" = "taxonomy7"))
abund_ubiq278_samp <- data.frame("SampleID"=c(colnames(abund_ubiq278_data[2:5]),
                                              colnames(abund_ubiq278_data[6:9])),
                                 "Food" = c("Honey","Honey","Honey","Honey",
                                            "Termites","Termites","Termites","Termites"),
                                 "Region" = c("ellioti","schweinfurthii","troglodytes","verus",
                                              "ellioti","schweinfurthii","troglodytes","verus"))
rownames(abund_ubiq278_taxa) <- abund_ubiq278_taxa$asvID
rownames(abund_ubiq278_samp) <- abund_ubiq278_samp$SampleID
tax.mat2 <- as.matrix(abund_ubiq278_taxa)
abund_ubiq278_data <- abund_ubiq278_data[,-1]
abund_ubiq278_data <- abund_ubiq278_data %>% mutate_if(is.character, as.numeric)
rownames(abund_ubiq278_data) <- abund_ubiq278$asvID
otu_tab2 <- otu_table(abund_ubiq278_data, taxa_are_rows = TRUE)
tax_tab2 <- tax_table(tax.mat2)
sam_tab2 <- sample_data(abund_ubiq278_samp)
ps2 <- phyloseq(otu_tab2, tax_tab2, sam_tab2, tree)
ps2 <- prune_taxa(no_mix, ps2)

circ2 <- ggtree(ps2, layout = "circular", branch.length = "none") + 
  geom_point(aes(x = x + hjust, color = Food, shape = Region)) +
  geom_cladelabel(node = 313, label = "Actinobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[1], hjust = "right", angle = 39) +
  geom_cladelabel(node = 333, label = "Bacteroidetes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[2], hjust = "left", angle = -58) +
  geom_cladelabel(node = 289, label = "Chloroflexi", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[3], hjust = "right", angle = -5) +
  geom_cladelabel(node = 285, label = "Cyanobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[4], hjust = "right", angle = -20) +
  geom_cladelabel(node = 423, label = "Fibrobacteres", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[5], hjust = "right", angle = 54) +
  geom_cladelabel(node = c(512, 448), label = "Firmicutes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[6], hjust = "left", angle = 79) +
  geom_cladelabel(node = 291, label = "Kiritimatiellaeota", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[7], hjust = "right", angle = 4) +
  geom_cladelabel(node = 277, label = "Proteobacteria", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5,  color = hue_pal()(9)[8], hjust = "right", angle = -14) +
  geom_cladelabel(node = 304, label = "Spirochaetes", fontsize = 3, offset = 3, offset.text = 0.2, barsize = 1.5, color = hue_pal()(9)[9], hjust = "right", angle = 21) +
  scale_color_manual(values = c("gold","red"),
                     na.translate = F) +
  scale_shape_manual(values = c(18, 16, 17, 15),
                     labels = c("N-C", "East", "Central", "West"),
                     na.translate = F) +
  labs(color = "Consumption of",
       shape = "In region") +
  geom_text(x = -45, y = -40, label = "b) Diet Indicator Taxa", size = 6) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme(legend.position = c(-0.07, 0.5),
        legend.background = element_blank(),
        plot.margin = unit(c(-0.5,-5,-0.5,-0.5), "cm")) +
  xlim(0,36)
circ2



# Multipanel
png(file = "PNGs/SupplementaryFigureS6.png", width = 14, height = 6, unit = "in", res = 300)
multiplot(circ1, circ2, ncol = 2)
dev.off()



############################## Main Analyses and Figures ######################################
################################### _1. Site Map ##############################################
# Basic site map
# The final map for Figure 1 with subspecies ranges and algae, honey, nuts, termites consumption was made in QGIS
sitemap <- dplyr::select(bac_euk$map_loaded,
                  subspecies, biom_AK, habitat, Site, long, lat)
sitemap <- sitemap[ !duplicated(sitemap$Site), ]
# For review, get plant richness by habitat type
plant_richness <- as.data.frame(specnumber(plantspcomp)) %>%
  rownames_to_column(var = "Site") %>%
  left_join(sitemap, by = "Site") %>%
  group_by(habitat) %>%
  summarise(mean = mean(`specnumber(plantspcomp)`),
            se = se(`specnumber(plantspcomp)`))

plant_richness <- as.data.frame(specnumber(plantspcomp)) %>%
  rownames_to_column(var = "Site") %>%
  left_join(biomAleman, by = c("Site" = "Site_Name2")) %>%
  group_by(biome_Aleman) %>%
  summarise(mean = mean(`specnumber(plantspcomp)`),
            se = se(`specnumber(plantspcomp)`))
plant_richness

# write.csv(sitemap, "sitecoords.csv")
sitediet <- select(bac_euk$map_loaded,
                  subspecies, biom_AK, Site, long, lat, algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', termites, tubers, water)
sitediet <- sitediet[ !duplicated(sitediet$Site), ]
# write.csv(sitediet, "sitediet.csv")
world <- map_data("world")
ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group),
               fill = "grey", color = "white", size = 0.25) +
  geom_point(data = sitemap, 
             aes(x = long, y = lat, color = subspecies, shape = biom_AK),
             size = 2) +
  geom_point(data = subset(sitemap, biom_AK == "forest"),
             aes(x = long, y = lat), size = 2, pch = 1, alpha = 0.75) +
  geom_point(data = subset(sitemap, biom_AK == "savanna"),
             aes(x = long, y = lat), size = 2, pch = 2, alpha = 0.75) +
  geom_segment(aes(x = 55, xend = 55, y = -35, yend = -31),
               arrow = arrow(length = unit(0.4, "cm")), colour = "black",
               inherit.aes = FALSE) +
  geom_text(aes(x = 55, y = -29, label = "N"), size = 5) +
  geom_segment(aes(x = 40.5, xend = 49.509, y = -35, yend = -35), colour = "black", size = 2) +
  geom_text(aes(x = 45, y = -33, label = "1000 km"), size = 3) +
  coord_fixed(1.3) +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Region",
       shape = "Habitat") +
  xlim(-20, 55) + 
  ylim(-35, 37) +
  theme(legend.position = c(0.22,0.25),
        legend.title = element_text(size = 10),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.spacing.y = unit(0.1, "cm"),
        panel.grid = element_blank())



################################### _2. Main Taxa #############################################
# a) Bacteria families by site/subsp b) Parasites by site/subsp
# a) Bacteria
tax_sum_families <- summarize_taxonomy(bac_euk, level = 5, report_higher_tax = FALSE)
fam <- as.data.frame(t(tax_sum_families))
plot_ts_heatmap(tax_sum_families, bac_euk$map_loaded, 0.01, "Site")
tax_sum_phyla <- summarize_taxonomy(bac_euk, level = 2, report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_phyla, bac_euk$map_loaded, 0.01, "Site") + coord_flip()
plot_ts_heatmap(tax_sum_phyla, bac_euk$map_loaded, 0.01, "Blastocystis_59") + coord_flip()
tax_sum_genera <- summarize_taxonomy(bac_euk, level = 6, report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_genera, bac_euk$map_loaded, 0.01, "Blastocystis_59") + coord_flip()
# With mean rel abund 0.01
fam <- select(fam,
              Atopobiaceae, Bacteroidales_RF16_group, Bifidobacteriaceae, Clostridiales_vadinBB60_group, Eggerthellaceae, Erysipelotrichaceae, Lachnospiraceae, Methanomethylophilaceae, Muribaculaceae, Prevotellaceae, Rikenellaceae, Ruminococcaceae, Spirochaetaceae, Veillonellaceae)
sort(colSums(fam))
fam$sampleID <- rownames(fam)
site_subspecies <- select(bac_euk$map_loaded, sampleID, Site, subspecies)
fam <- inner_join(fam, site_subspecies, by = c("sampleID" = "sampleID"))
fam_long <- melt(fam, measure.vars = c("Atopobiaceae",
                                       "Bacteroidales_RF16_group", 
                                       "Bifidobacteriaceae", 
                                       "Clostridiales_vadinBB60_group", 
                                       "Eggerthellaceae", 
                                       "Erysipelotrichaceae", 
                                       "Lachnospiraceae", 
                                       "Methanomethylophilaceae", 
                                       "Muribaculaceae", 
                                       "Prevotellaceae", 
                                       "Rikenellaceae", 
                                       "Ruminococcaceae", 
                                       "Spirochaetaceae", 
                                       "Veillonellaceae"))
site_abund <- ddply(fam_long, c("subspecies","Site","variable"), summarise,
                    abund = round((mean(value)*100), digits = 1))
site_abund <- mutate_(site_abund, scaled = ~ scales::rescale(abund))

# Order sites
site_abund$Site <- factor(site_abund$Site, 
                          levels = c("Mt_Cameroon","Gashaka","Gishwati","Nyungwe","Issa","Budongo","Bwindi","Conkouati","Goualougo","Loango","Lope","Mts_de_Cristal","Bakoun","Sangaredi","Sobeya","Sobory","Boe","Comoe_GEPRENAF","Djouroutou","Mt_Sangbe","Tai_Ecotourism","Tai_Recherche","East_Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi"))

# Update for West to east
site_abund$Site <- factor(site_abund$Site, 
                          levels = c("Bakoun","Sangaredi","Sobeya","Sobory","Boe","Comoe_GEPRENAF","Djouroutou","Mt_Sangbe","Tai_Ecotourism","Tai_Recherche","East_Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi","Mt_Cameroon","Gashaka","Conkouati","Goualougo","Loango","Lope","Mts_de_Cristal","Gishwati","Nyungwe","Issa","Budongo","Bwindi"))

# Update for clustering?
bac_clust_data <- plot_taxa_bars(tax_sum_families,
                                 bac_euk$map_loaded,
                                 "Site",
                                 14,
                                 data_only = TRUE) %>%
  dcast(group_by ~ taxon) %>%
  column_to_rownames(var = "group_by") %>%
  select(-`NA`, -Other)
bac_bc_fam <- vegdist(bac_clust_data, method = "bray")
bac_clust <- hclust(bac_bc_fam, method = "ward.D2")
plot(bac_clust) # Doesn't group by region
library(dendextend)
bc_den <- as.matrix(bacteria.distance)
rownames(bc_den) <- bac_euk$map_loaded$subspecies
colnames(bc_den) <- bac_euk$map_loaded$subspecies
bc_den <- as.dist(bc_den)
bac_clust <- hclust(bc_den, method = "ward.D2")
den <- as.dendrogram(bac_clust)
colors_to_use <- as.numeric(as.factor(bac_euk$map_loaded$subspecies))
colors_to_use <- colors_to_use[order.dendrogram(den)]
labels_colors(den) <- colors_to_use
labels_cex(den) <- 0.1
par(oma = c(0,0,0,0),
    mar = c(2,5,0,0))
plot_horiz.dendrogram(den) # Mostly grouped by region but not perfect. Let's argue not to do it this way. We can still sort the y axis though according to the subspecies effect.

# Subspecies effect
kw <- taxa_summary_by_sample_type(tax_sum_families, bac_euk$map_loaded, 
                                  type_header = 'subspecies', 
                                  filter_level = 0.01, 
                                  test_type = 'KW') %>%
  rownames_to_column(var = "fam") %>%
  mutate(pvalsBonRd = round(pvalsBon, digits = 4)) %>%
  filter(fam == "Atopobiaceae" |
           fam == "Bacteroidales_RF16_group" |
           fam == "Bifidobacteriaceae" |
           fam == "Clostridiales_vadinBB60_group" |
           fam == "Eggerthellaceae" |
           fam == "Erysipelotrichaceae" |
           fam == "Lachnospiraceae" |
           fam == "Methanomethylophilaceae" |
           fam == "Muribaculaceae" |
           fam == "Prevotellaceae" | 
           fam == "Rikenellaceae" |
           fam == "Ruminococcaceae" |
           fam == "Spirochaetaceae" |
           fam == "Veillonellaceae")
kw$fam <- rownames(kw)
kw$pvalsBonRd <- round(kw$pvalsBon, digits = 4)

# Top family mean ± SE
tax_sum_families_t <- as.data.frame(t(tax_sum_families))
mean(tax_sum_families_t$Prevotellaceae)
se(tax_sum_families_t$Prevotellaceae)
mean(tax_sum_families_t$Lachnospiraceae)
se(tax_sum_families_t$Lachnospiraceae)
mean(tax_sum_families_t$Ruminococcaceae)
se(tax_sum_families_t$Ruminococcaceae)
mean(tax_sum_families_t$Erysipelotrichaceae)
se(tax_sum_families_t$Erysipelotrichaceae)
mean(tax_sum_families_t$Rikenellaceae)
se(tax_sum_families_t$Rikenellaceae)
top_families <- site_abund %>%
  group_by(variable) %>%
  summarise(mean = mean(abund),
            se = se(abund)) %>%
  arrange(mean) %>%
  mutate(variable = as.character(variable))
site_abund$variable_sorted <- factor(site_abund$variable,
                                     levels = top_families$variable)

pdf(file = "forPPT3bac.pdf", width = 7, height = 3)
ggplot(site_abund, aes(x = Site, y = variable_sorted, fill = scaled)) +
  geom_tile(color = "black", size = 0.25) +
  geom_text(data = site_abund, aes_string(label = "abund"), size = 1.75) +
  scale_fill_gradientn(colours = c('blue', 'white', 'red'),
                       values = c(
                         min(site_abund$scaled),
                         mean(site_abund$scaled),
                         max(site_abund$scaled)
                       )) +
  scale_x_discrete(expand = c(0,0),
                   labels = c("Bakoun","Sangaredi","Sobeya","Sobory","Boé","Comoe-GEPRENAF","Djouroutou","Mt. Sangbé","Taï Ecotourism","Taï Recherche","East Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi","Mt. Cameroon","Gashaka","Conkouati","Goualougo","Loango","Lopé","Mts. de Cristal","Gishwati","Nyungwe","Issa","Budongo","Bwindi")) +
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(c(0,-1.5,0,0))),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, 
                                   margin = margin(c(-1.5,0,0,0))))
dev.off()


# b) Parasites
# Heatmaps for prevalence, by Site (ordered by subspecies)
site_prev_num <- melt(bac_euk$map_loaded, measure.vars = c("Blastocystis_59",
                                                           "Blepharocorys_uncinata",
                                                           "Chilomastix_mesnili",
                                                           "Dientamoeba_fragilis",
                                                           "Entamoeba_hartmanni",
                                                           "Entamoeba_muris",
                                                           "Haemonchus_contortus",
                                                           "Strongyloides_fuelleborni",
                                                           "Tetratrichomonas_2",
                                                           "Tetratrichomonas_31",
                                                           "Trichomonadidae_12_51",
                                                           "Trichomonadidae_15",
                                                           "Trichomonadidae_8",
                                                           "Troglodytella_abrassarti"))
site_prev <- ddply(site_prev_num, c("Site","variable"), summarise,
                   prevalence = round((sum(value)/length(value))*100, digits = 0)) %>%
  mutate_(scaled = ~ scales::rescale(prevalence))

# Update for West to east
site_prev$Site <- factor(site_prev$Site, 
                          levels = c("Bakoun","Sangaredi","Sobeya","Sobory","Boe","Comoe_GEPRENAF","Djouroutou","Mt_Sangbe","Tai_Ecotourism","Tai_Recherche","East_Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi","Mt_Cameroon","Gashaka","Conkouati","Goualougo","Loango","Lope","Mts_de_Cristal","Gishwati","Nyungwe","Issa","Budongo","Bwindi"))

# Update for ordering rows
top_parasites <- bac_euk$map_loaded %>%
  select(Blastocystis_59,
         Blepharocorys_uncinata,
         Chilomastix_mesnili,
         Dientamoeba_fragilis,
         Entamoeba_hartmanni,
         Entamoeba_muris,
         Haemonchus_contortus,
         Strongyloides_fuelleborni,
         Tetratrichomonas_2,
         Tetratrichomonas_31,
         Trichomonadidae_12_51,
         Trichomonadidae_15,
         Trichomonadidae_8,
         Troglodytella_abrassarti) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parasite") %>%
  mutate(NumPres = rowSums(.[,2:561])) %>%
  arrange(NumPres) %>%
  select(Parasite, NumPres)
site_prev$variable_sorted <- factor(site_prev$variable,
                                    levels = top_parasites$Parasite)

pdf(file = "forPPT3euk.pdf", width = 7, height = 3)
ggplot(site_prev, aes(x = Site, y = variable_sorted, fill = scaled)) +
  geom_tile(color = "black", size = 0.25) +
  geom_text(data = site_prev, aes_string(label = "prevalence"), size = 1.75) +
  scale_fill_gradientn(colours = c('blue', 'white', 'red'),
                       values = c(
                         min(site_prev$scaled),
                         mean(site_prev$scaled),
                         max(site_prev$scaled)
                       )) +
  scale_x_discrete(expand = c(0,0),
                   labels = c("Bakoun","Sangaredi","Sobeya","Sobory","Boé","Comoe-GEPRENAF","Djouroutou","Mt. Sangbé","Taï Ecotourism","Taï Recherche","East Nimba","Grebo","Sapo","Bafing","Dindefelo","Kayan","Outamba-Kilimi","Mt. Cameroon","Gashaka","Conkouati","Goualougo","Loango","Lopé","Mts. de Cristal","Gishwati","Nyungwe","Issa","Budongo","Bwindi")) +
  scale_y_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10, face = "italic", margin = margin(c(0,0,0,0))),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, 
                                   margin = margin(c(-1.5,0,0,0))))
dev.off()

# Logistic regressions for region effect
Anova(glm(Blastocystis_59 ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Blepharocorys_uncinata ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Chilomastix_mesnili ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Dientamoeba_fragilis ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Entamoeba_hartmanni ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Entamoeba_muris ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Haemonchus_contortus ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Strongyloides_fuelleborni ~ subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Tetratrichomonas_2~subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Tetratrichomonas_31~subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Trichomonadidae_12_51~subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Trichomonadidae_15~subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Trichomonadidae_8~subspecies, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Troglodytella_abrassarti ~ subspecies, family = binomial, data = bac_euk$map_loaded))

# Logistic regressions for region, sex, human footprint effect
Anova(glm(Blastocystis_59 ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Blepharocorys_uncinata ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Chilomastix_mesnili ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Dientamoeba_fragilis ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Entamoeba_hartmanni ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Entamoeba_muris ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Haemonchus_contortus ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Strongyloides_fuelleborni ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Tetratrichomonas_2~subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Tetratrichomonas_31~subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Trichomonadidae_12_51~subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Trichomonadidae_15~subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glm(Trichomonadidae_8~subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))
Anova(glm(Troglodytella_abrassarti ~ subspecies + Sex + Human_Footprint, family = binomial, data = bac_euk$map_loaded))

# Mixed Effects Logistic regressions for region, sex, human footprint effect
Anova(glmer(Blastocystis_59 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Blepharocorys_uncinata ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Chilomastix_mesnili ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glmer(Dientamoeba_fragilis ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # No longer sig
Anova(glmer(Entamoeba_hartmanni ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Entamoeba_muris ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Haemonchus_contortus ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # No longer sig
Anova(glmer(Strongyloides_fuelleborni ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Tetratrichomonas_2 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Tetratrichomonas_31 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glmer(Trichomonadidae_12_51 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glmer(Trichomonadidae_15 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded)) # Sig
Anova(glmer(Trichomonadidae_8 ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))
Anova(glmer(Troglodytella_abrassarti ~ habitat + Sex + Human_Footprint + (1|subspecies), family = binomial, data = bac_euk$map_loaded))

# Test prevalence and climate
site_prev_clim <- inner_join(site_prev, site_info, by = "Site")
summary(lm(site_prev_clim$prevalence ~ site_prev_clim$AnPercip))
ggplot(site_prev_clim, aes(AnPercip, prevalence, colour = variable)) +
  geom_point() +
  geom_smooth(method = lm)
# 3,4,7,10,11,12 glm
# 3, 10, 11, 12 glmer
ggplot(site_prev_clim, aes(Human_Footprint, prevalence, colour = variable)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "Human Footprint",
       y = "Site Prevalence",
       colour = "Parasite") +
  geom_smooth(data = subset(site_prev_clim, 
                              variable == "Chilomastix_mesnili" |
                              variable == "Tetratrichomonas_31" |
                              variable == "Trichomonadidae_12_51" |
                              variable == "Trichomonadidae_15"),
                            method = lm, se = FALSE)
sf <- subset(site_prev_clim, variable == "Strongyloides_fuelleborni")
summary(lm(sf$prevalence ~ sf$AnPercip))



################################### _3. PCoAs #################################################
# Multipanel- a) Genetics/subspecies PCoA, b) Bacteria/subsp. PCoA c) Parasite/subsp. PCoA
# a) Genetics/subspecies PCoA
gen_mat_df <- select(gen_mat_df, -site)
genetic.distance <- as.dist(gen_mat_df)
pcoa.gen <- cmdscale(genetic.distance, k = nrow(gen_mat_df)-1, eig=T)
eigenvals(pcoa.gen)/sum(eigenvals(pcoa.gen)) # 53.3, 17.5
pcoa_gen_df <- data.frame("Axis01" = scores(pcoa.gen)[,1],
                          "Axis02" = scores(pcoa.gen)[,2])
pcoa_gen_df$Site <- rownames(pcoa_gen_df)
pcoa_gen_df <- left_join(pcoa_gen_df, sitemap,
                         by = c("Site" = "Site"))
# Add subspecies info from Jack paper even if site not in our sampling. 
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Bakoun-Sobory" |
      pcoa_gen_df$Site[i] == "Comoe_East") {
    pcoa_gen_df$subspecies[i] <- "verus"
  }
}
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Korup") {
    pcoa_gen_df$subspecies[i] <- "ellioti"
  }
}
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "La_Belgique") {
    pcoa_gen_df$subspecies[i] <- "troglodytes"
  }
}
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Rubi_Tele" |
      pcoa_gen_df$Site[i] == "Bili" |
      pcoa_gen_df$Site[i] == "Chinko" |
      pcoa_gen_df$Site[i] == "Ngogo" |
      pcoa_gen_df$Site[i] == "Issa") {
    pcoa_gen_df$subspecies[i] <- "schweinfurthii"
  }
}
# Add biom_AK for genetic sites
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Korup" |
      pcoa_gen_df$Site[i] == "La_Belgique" |
      pcoa_gen_df$Site[i] == "Rubi_Tele" |
      pcoa_gen_df$Site[i] == "Bili" |
      pcoa_gen_df$Site[i] == "Chinko"|
      pcoa_gen_df$Site[i] == "Ngogo") {
    pcoa_gen_df$biom_AK[i] <- "forest"
  }
}
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Bakoun-Sobory" |
      pcoa_gen_df$Site[i] == "Comoe_East") {
    pcoa_gen_df$biom_AK[i] <- "savanna"
  }
}
# Add habitat for genetic sites
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Korup" |
      pcoa_gen_df$Site[i] == "La_Belgique" |
      pcoa_gen_df$Site[i] == "Rubi_Tele" |
      pcoa_gen_df$Site[i] == "Bili" |
      pcoa_gen_df$Site[i] == "Chinko"|
      pcoa_gen_df$Site[i] == "Ngogo") {
    pcoa_gen_df$habitat[i] <- "forest"
  }
}
for (i in 1:nrow(pcoa_gen_df)) {
  if (pcoa_gen_df$Site[i] == "Bakoun-Sobory" |
      pcoa_gen_df$Site[i] == "Comoe_East") {
    pcoa_gen_df$habitat[i] <- "savanna"
  }
}

# Add new Aleman habitat
pcoa_gen_df2 <- left_join(pcoa_gen_df, biomAleman, by = c("Site" = "Site_Name2"))
micro.hulls <- ddply(pcoa_gen_df2, "subspecies", find_hull)
ggplot(pcoa_gen_df, aes(Axis01, Axis02, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1: 53.3% Variation Explained") +
  ylab("PC2: 17.5% Variation Explained") +
  labs(colour = "Subspecies") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))

# b) Bacteria/subspecies PCoA
pcoa.bac <- cmdscale(bacteria.distance, k = nrow(bac_euk$map_loaded)-1, eig=T)
eigenvals(pcoa.bac)/sum(eigenvals(pcoa.bac)) # 18.9, 6.5
bac_euk$map_loaded$Axis01 <- scores(pcoa.bac)[,1]
bac_euk$map_loaded$Axis02 <- scores(pcoa.bac)[,2]
micro.hulls <- ddply(bac_euk$map_loaded, "subspecies", find_hull)
ggplot(bac_euk$map_loaded, aes(Axis01, Axis02, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1: 18.9% Variation Explained") +
  ylab("PC2: 6.5% Variation Explained") +
  labs(colour = "Subspecies") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))

# c) Parasite/subspecies PCoA
pcoa.par <- cmdscale(parasites.distance, k = nrow(bac_euk$map_loaded)-1, eig=T)
eigenvals(pcoa.par)/sum(eigenvals(pcoa.par)) # 25.8, 17.8
bac_euk$map_loaded$Axis01p <- scores(pcoa.par)[,1]
bac_euk$map_loaded$Axis02p <- scores(pcoa.par)[,2]
micro.hulls <- ddply(bac_euk$map_loaded, "subspecies", find_hullp)
ggplot(bac_euk$map_loaded, aes(Axis01p, Axis02p, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1: 32.1% Variation Explained") +
  ylab("PC2: 18.5% Variation Explained") +
  labs(colour = "Subspecies") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))

# Make dataframe for 3 panel
pcoa_gen_df$dataset <- "a"
pcoa_multi_a <- select(pcoa_gen_df, Axis01, Axis02, subspecies, dataset, Site)
pcoa_multi_b <- select(bac_euk$map_loaded, Axis01, Axis02, subspecies, Site)
pcoa_multi_b$dataset <- "b"
pcoa_multi_c <- select(bac_euk$map_loaded, Axis01p, Axis02p, subspecies)
pcoa_multi_c$dataset <- "c"
names(pcoa_multi_c) <- c("Axis01","Axis02","subspecies","dataset")
pcoa_multi_long <- rbind(pcoa_multi_a,pcoa_multi_b,pcoa_multi_c)
pcoa_multi_long$subspecies <- factor(pcoa_multi_long$subspecies,
                                     levels = c("verus", "ellioti", 
                                                "troglodytes","schweinfurthii"))
pcoa_multi_long$dataset <- as.factor(pcoa_multi_long$dataset)
micro.hulls <- ddply(pcoa_multi_long, c("dataset","subspecies"), find_hull)
facet_names <- c('a' = "a) Genetic",
                 'b' = "b) Bacteria/Archaea",
                 'c' = "c) Parasites")
label.df.2 <- data.frame(dataset = c("a","b","c"),
                         x = c(0.28, 0.27, 0.27),
                         y = c(-0.24, -0.33, -0.48),
                         label = c("53.3%, 17.5%", "18.9%, 6.5%", "25.8%, 17.8%"))
ggplot(pcoa_multi_long, aes(Axis01, Axis02, colour = subspecies, fill = subspecies)) +
  geom_polygon(data = micro.hulls, alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.4) +
  geom_text(data = label.df.2, aes(x = x, y = y, label = label), size = 3, 
            inherit.aes = F) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00")) +
  scale_fill_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00")) +
  facet_wrap(~dataset, ncol = 3, scales = "free", labeller = as_labeller(facet_names)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.margin = margin(t = -0.3, unit = "cm"),
        legend.key.size = unit(0.25, "cm"),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(c(0,-2,0,0))), 
        axis.title.x = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))

### Update - Also do PCoA on climate and distance!
sitemap$sampleID <- rownames(sitemap)
# d) Climate (Update, this is going to supplement)
pcoa.cli <- cmdscale(sl.climate.distance, k = 28, eig=T)
eigenvals(pcoa.cli)/sum(eigenvals(pcoa.cli)) # 62.4, 25.8
pcoa_clim_df <- data.frame("Axis01" = scores(pcoa.cli)[,1],
                           "Axis02" = scores(pcoa.cli)[,2])
pcoa_clim_df$sampleID <- rownames(pcoa_clim_df)
pcoa_clim_df <- left_join(pcoa_clim_df, sitemap,
                          by = c("sampleID" = "sampleID"))
micro.hulls <- ddply(pcoa_clim_df, "subspecies", find_hull)
pdf(file = "forPPTclim1.pdf", width = 6, height = 4)
ggplot(pcoa_clim_df, aes(Axis01, Axis02, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = "PC1: 62.4% Variation Explained",
       y = "PC2: 25.8% Variation Explained",
       colour = "Region") +
  scale_color_discrete(breaks = c("verus", "ellioti", "troglodytes", "schweinfurthii"),
                       labels = c("West", "Nigeria-Cameroon", "Central", "East")) +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()

# e) Diet
pcoa.die <- cmdscale(sl.diet.distance, k = 28, eig=T)
eigenvals(pcoa.die)/sum(eigenvals(pcoa.die)) # 27.7, 23.1
pcoa_diet_df <- data.frame("Axis01" = scores(pcoa.die)[,1],
                          "Axis02" = scores(pcoa.die)[,2])
pcoa_diet_df$sampleID <- rownames(pcoa_diet_df)
pcoa_diet_df <- left_join(pcoa_diet_df, sitemap,
                         by = c("sampleID" = "sampleID"))
micro.hulls <- ddply(pcoa_diet_df, "subspecies", find_hull)
ggplot(pcoa_diet_df, aes(Axis01, Axis02, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1: 29.9% Variation Explained") +
  ylab("PC2: 24.3% Variation Explained") +
  labs(colour = "Subspecies") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))

# New 5 panel graph
pcoa_mult_a <- select(pcoa_gen_df, Axis01, Axis02, subspecies, biom_AK, habitat, Site)
pcoa_mult_a$dataset <- "a"
pcoa_mult_b <- select(pcoa_clim_df, Axis01, Axis02, subspecies, biom_AK, habitat, Site)
pcoa_mult_b$dataset <- "b"
pcoa_mult_c <- select(pcoa_diet_df, Axis01, Axis02, subspecies, biom_AK, habitat, Site)
pcoa_mult_c$dataset <- "c"
pcoa_mult_d <- select(bac_euk$map_loaded, Axis01, Axis02, subspecies, biom_AK, habitat, Site)
pcoa_mult_d$dataset <- "d"
pcoa_mult_e <- select(bac_euk$map_loaded, Axis01p, Axis02p, subspecies, biom_AK, habitat, Site)
pcoa_mult_e$dataset <- "e"
pcoa_mult_e$Axis01p <- -pcoa_mult_e$Axis01p # Flip Axis 1
pcoa_mult_e$Axis02p <- -pcoa_mult_e$Axis02p # Flip Axis 2
names(pcoa_mult_e) <- c("Axis01","Axis02","subspecies","biom_AK","habitat","Site","dataset")
pcoa_mult5_long <- rbind(pcoa_mult_a,pcoa_mult_b,pcoa_mult_c,pcoa_mult_d,pcoa_mult_e)
pcoa_mult5_long$subspecies <- factor(pcoa_mult5_long$subspecies,
                                    levels = c("verus","ellioti","troglodytes","schweinfurthii"))
pcoa_mult5_long$dataset <- as.factor(pcoa_mult5_long$dataset)
micro.hulls <- ddply(pcoa_mult5_long, c("dataset","subspecies"), find_hull)
facet_names.pcoa <- c('a' = "a) Genetic",
                      'b' = "b) Climate",
                      'c' = "c) Diet/Tool Use",
                      'd' = "d) Bacteria/Archaea",
                      'e' = "e) Parasites")
label.df.pcoa <- data.frame(dataset = c("a","b","c","d","e"),
                            x = c(-0.14,-1.8,-0.5, 0.25, 0.23),
                            y = c(-0.24, -1.75, -0.6, -0.33, -0.48),
                            label = c("53.3%,17.5%","62.4%,25.8%","29.9%,24.3%",
                                      "18.9%,6.5%","25.8%, 17.8%"))
ggplot(pcoa_mult5_long, 
       aes(Axis01, Axis02, colour = subspecies, fill = subspecies)) +
  geom_polygon(data = micro.hulls, alpha = 0.1, show.legend = F) +
  geom_point(size = 1, alpha = 0.4, aes(shape = habitat)) +
  geom_text(data = label.df.pcoa, aes(x = x, y = y, label = label), size = 3, 
            inherit.aes = F) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00")) +
  scale_fill_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00")) +
  scale_shape_manual(values = c(24, 25, 21),
                     labels = c("forest", "forest (low density)", "savanna")) +
  facet_wrap(~dataset, ncol = 3, scales = "free", labeller = as_labeller(facet_names.pcoa)) +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)),
         shape = guide_legend(override.aes = list(alpha=1, size = 2))) +
  theme(legend.position = c(0.85, 0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(0, "cm"),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(c(0,-2,0,0))), 
        axis.title.x = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))



#### _Update: Vegetation instead of climate, no diet ####
rn2 <- rownames(plantmeta)
plantmeta <- left_join(plantmeta, site_info, by = c("Site" = "Site"))
rownames(plantmeta) <- rn2
eigenvals(pcoa.veg)/sum(eigenvals(pcoa.veg)) # 21.1, 11.4
pcoa_veg_df <- data.frame("Axis01" = scores(pcoa.veg)[,1],
                          "Axis02" = scores(pcoa.veg)[,2])
pcoa_veg_df$Site <- rownames(pcoa_veg_df)
pcoa_veg_df <- left_join(pcoa_veg_df, sitemap,
                         by = c("Site" = "Site"))
micro.hulls <- ddply(pcoa_veg_df, "subspecies", find_hull)
ggplot(pcoa_veg_df, aes(Axis01, Axis02, colour = subspecies)) +
  geom_polygon(data = micro.hulls, aes(colour = subspecies, fill = subspecies),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  xlab("PC1: 21.1% Variation Explained") +
  ylab("PC2: 11.4% Variation Explained") +
  labs(colour = "Subspecies") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14))

pcoa_mult_b <- select(pcoa_veg_df, Axis01, Axis02, subspecies, biom_AK, habitat, Site)
pcoa_mult_b$dataset <- "b"

# Update to remove diet
pcoa_mult5_long <- rbind(pcoa_mult_a,pcoa_mult_b,pcoa_mult_d,pcoa_mult_e)
pcoa_mult5_long$subspecies <- factor(pcoa_mult5_long$subspecies,
                                     levels = c("verus","ellioti","troglodytes","schweinfurthii"))
pcoa_mult5_long$dataset <- as.factor(pcoa_mult5_long$dataset)
pcoa_mult5_long <- left_join(pcoa_mult5_long, biomAleman, by = c("Site" = "Site_Name2"))
pcoa_mult5_long$biome_Aleman <- factor(pcoa_mult5_long$biome_Aleman,
                                       levels = c("forest", "forest mosaic",
                                                  "savanna mosaic", "savanna"))
micro.hulls <- ddply(pcoa_mult5_long, c("dataset","subspecies"), find_hull)
facet_names.pcoa <- c('a' = "a) Genetic",
                      'b' = "b) Vegetation",
                      'd' = "c) Prokaryotes",
                      'e' = "d) Parasites")
label.df.pcoa <- data.frame(dataset = c("a","b","d","e"),
                            x = c(-0.14,-0.24, 0.25, 0.23),
                            y = c(-0.23, -0.59, -0.32, -0.47),
                            label = c("53.3%,17.5%","62.4%,25.8%",
                                      "18.9%,6.5%","25.8%,17.8%"))
# pdf(file = "Figure3.pdf", width = 7, height = 4)
tiff(file = "PNGs/Figure3.tiff", width = 6.87, height = 4, units = "in", res = 300)
ggplot(pcoa_mult5_long, 
       aes(Axis01, Axis02, colour = subspecies, fill = subspecies)) +
  geom_polygon(data = micro.hulls, alpha = 0.1, show.legend = F) +
  geom_point(size = 1, alpha = 0.4, aes(shape = biome_Aleman)) +
  geom_text(data = label.df.pcoa, aes(x = x, y = y, label = label), size = 3, 
            inherit.aes = F) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00"),
                     labels = c("West", "Nigeria-Cameroon", "Central", "East")) +
  scale_fill_manual(values = c("#C77CFF","#F8766D","#00BFC4","#7CAE00"),
                    guide = FALSE) +
  scale_shape_manual(values = c(24, 25, 21, 22),
                     labels = c("forest", "forest mosaic", "savanna mosaic", "savanna")) +
  facet_wrap(~dataset, ncol = 2, scales = "free", labeller = as_labeller(facet_names.pcoa)) +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)),
         shape = guide_legend(override.aes = list(alpha=1, size = 2))) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.spacing.y = unit(0, "cm"),
        legend.box.margin = margin(0,0,0,-15),
        legend.background = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(c(0,-2,0,0))), 
        axis.title.x = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))
dev.off()



############################### _4. 16S/18S Mantel Tests ######################################
# Run section 5 first to make some dataframes!
# Figure 4
mantel(parasites.distance, bacteria.distance) # p = 0.001, r = 0.55
mantel(elli.parasites.distance, elli.bacteria.distance) # p = 0.001, r = 0.81
mantel(schwein.parasites.distance, schwein.bacteria.distance) # p = 0.001, r = 0.24
mantel(trog.parasites.distance, trog.bacteria.distance) # p = 0.001, r = 0.48
mantel(verus.parasites.distance, verus.bacteria.distance) # p = 0.001, r = 0.36
# Make 10 dataframes, combine, make 5 column facet wrap
# Already have 8, just make elli
elli.bacteria.distance.mat <- as.matrix(elli.bacteria.distance)
elli.bacteria.distance.mat[upper.tri(elli.bacteria.distance.mat, diag = TRUE)] <- NA
elli.bacteria.distance.df <- as.data.frame(elli.bacteria.distance.mat)
elli.bacteria.distance.df$sampleID <- rownames(elli.bacteria.distance.df)
elli.bacteria.distance.df.long <- melt(elli.bacteria.distance.df, id.vars = "sampleID")
elli.bacteria.distance.df.long <- na.omit(elli.bacteria.distance.df.long)
names(elli.bacteria.distance.df.long) <- c("sampleID","variable","BC")
elli.bacteria.distance.df.long$subspecies <- "ellioti"
elli.parasites.distance.mat <- as.matrix(elli.parasites.distance)
elli.parasites.distance.mat[upper.tri(elli.parasites.distance.mat, diag = TRUE)] <- NA
elli.parasites.distance.df <- as.data.frame(elli.parasites.distance.mat)
elli.parasites.distance.df$sampleID <- rownames(elli.parasites.distance.df)
elli.parasites.distance.df.long <- melt(elli.parasites.distance.df, id.vars = "sampleID")
elli.parasites.distance.df.long <- na.omit(elli.parasites.distance.df.long)
names(elli.parasites.distance.df.long) <- c("sampleID","variable","Jac")
elli.parasites.distance.df.long$subspecies <- "ellioti"

bcj1 <- bacteria.distance.df.long[,1:3]
names(bcj1) <- c("sampleID","variable","bray")
bcj1$jac <- parasites.distance.df.long$Jac
bcj1$dataset <- "all"
bcj2 <- elli.bacteria.distance.df.long[,1:3]
names(bcj2) <- c("sampleID","variable","bray")
bcj2$jac <- elli.parasites.distance.df.long$value
bcj2$dataset <- "ellioti"
bcj3 <- schwein.bacteria.distance.df.long[,1:3]
names(bcj3) <- c("sampleID","variable","bray")
bcj3$jac <- schwein.parasites.distance.df.long$Jac
bcj3$dataset <- "schweinfurthii"
bcj4 <- trog.bacteria.distance.df.long[,1:3]
names(bcj4) <- c("sampleID","variable","bray")
bcj4$jac <- trog.parasites.distance.df.long$Jac
bcj4$dataset <- "troglodytes"
bcj5 <- verus.bacteria.distance.df.long[,1:3]
names(bcj5) <- c("sampleID","variable","bray")
bcj5$jac <- verus.parasites.distance.df.long$Jac
bcj5$dataset <- "verus"
bcj <- rbind(bcj1, bcj2, bcj3, bcj4, bcj5)
bcj$dataset <- as.factor(bcj$dataset)

# Graph
facet.names.bcj <- c('all' = "a) all",
                     'ellioti' = "c) N-C",
                     'schweinfurthii' = "e) East",
                     'troglodytes' = "d) Central",
                     'verus' = "b) West")
label.df.bcj <- data.frame(subspecies = c("all","ellioti","schweinfurthii","troglodytes","verus"),
                           dataset = c("all","ellioti","schweinfurthii","troglodytes","verus"),
                           x = c(0.5, 0.5, 0.5, 0.5, 0.5),
                           y = c(0.125,0.125,0.125,0.125,0.125),
                           label = c("n = 560, r = 0.55, p = 0.001",
                                     "n = 28, r = 0.81, p = 0.001",
                                     "n = 134, r = 0.24, p = 0.001",
                                     "n = 86, r = 0.48, p = 0.001",
                                     "n = 312, r = 0.36, p = 0.001"))
levels(bcj$dataset)
bcj$dataset <- factor(bcj$dataset,
                      levels = c("all", "verus", "ellioti", "troglodytes", "schweinfurthii"))
# pdf(file = "Figure4.pdf", width = 7.5, height = 3)
png(file = "PNGs/Figure4.png", width = 7.5, height = 3, units = "in", res = 300)
ggplot(data = bcj, aes(jac, bray)) +
  geom_jitter(size = 0.5, alpha = 0.05) +
  geom_smooth(se = FALSE) +
  geom_text(data = label.df.bcj, aes(x = x, y = y, label = label), size = 2.5) +
  labs(x = "Parasite Jaccard Dissimilarity",
       y = "Prokaryote Bray-Curtis Dissimilarity") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0,1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ dataset, ncol = 5, labeller = as_labeller(facet.names.bcj)) +
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 11))
dev.off()



##################### _5. Bray-Curtis: Geography, Climate, Diet ###############################
##################### __Part 1: Chimps ########################################################
# 12 panels (3 variables * all, schwein, trog, verus)
# No elli because only 2 sites (Geo, Climate, Diet the same), same as within vs among sites
# Not done at site level, done on whole dataset!
# Mantel Tests
mantel(bacteria.distance, geography.distance) # p = 0.001, r = 0.68
mantel(elli.bacteria.distance, elli.geography.distance) # p = 0.001, r = 0.66 (same for all)
mantel(schwein.bacteria.distance, schwein.geography.distance) # p = 0.001, r = 0.51
mantel(trog.bacteria.distance, trog.geography.distance) # p = 0.001, r = 0.52
mantel(verus.bacteria.distance, verus.geography.distance) # p = 0.001, r = 0.29

mantel(bacteria.distance, climate.distance) # p = 0.001, r = 0.40
mantel(elli.bacteria.distance, elli.climate.distance) # 0.001, 0.66
mantel(schwein.bacteria.distance, schwein.climate.distance) # p = 0.001, r = 0.58
mantel(trog.bacteria.distance, trog.climate.distance) # p = 0.001, r = 0.53
mantel(verus.bacteria.distance, verus.climate.distance) # p = 0.001, r = 0.24
mantel.partial(bacteria.distance, climate.distance, geography.distance) # p = 0.03, r = 0.03
mantel.partial(elli.bacteria.distance, elli.climate.distance, elli.geography.distance) # NA
mantel.partial(schwein.bacteria.distance, 
               schwein.climate.distance, 
               schwein.geography.distance) # p = 0.001, r = 0.34
mantel.partial(trog.bacteria.distance, 
               trog.climate.distance, 
               trog.geography.distance) # p = 0.003, r = 0.20
mantel.partial(verus.bacteria.distance, 
               verus.climate.distance, 
               verus.geography.distance) # p = 0.22, r = 0.02

mantel(bacteria.distance, diet.distance) # p = 0.001, r = 0.33
mantel(elli.bacteria.distance, elli.diet.distance) # p = 0.001, r = 0.66
mantel(schwein.bacteria.distance, schwein.diet.distance) # p = 0.001, r = 0.60
mantel(trog.bacteria.distance, trog.diet.distance) # p = 0.001, r = 0.30
mantel(verus.bacteria.distance, verus.diet.distance) # p = 0.001, r = 0.29
mantel.partial(bacteria.distance, diet.distance, geography.distance) # p = 0.001, r = 0.19
mantel.partial(elli.bacteria.distance, elli.diet.distance, elli.geography.distance) # NA
mantel.partial(schwein.bacteria.distance, 
               schwein.diet.distance, 
               schwein.geography.distance) # p = 0.001, r = 0.37
mantel.partial(trog.bacteria.distance, 
               trog.diet.distance, 
               trog.geography.distance) # p = 0.15, r = 0.06
mantel.partial(verus.bacteria.distance, 
               verus.diet.distance, 
               verus.geography.distance) # p = 0.001, r = 0.13

# Dataframes
# All
bacteria.distance.mat <- as.matrix(bacteria.distance)
bacteria.distance.mat[upper.tri(bacteria.distance.mat, diag = TRUE)] <- NA
bacteria.distance.df <- as.data.frame(bacteria.distance.mat)
bacteria.distance.df$sampleID <- rownames(bacteria.distance.df)
bacteria.distance.df.long <- melt(bacteria.distance.df, id.vars = "sampleID")
bacteria.distance.df.long <- na.omit(bacteria.distance.df.long)
geography.distance.mat[upper.tri(geography.distance.mat, diag = TRUE)] <- NA
geography.distance.df <- as.data.frame(geography.distance.mat)
geography.distance.df$sampleID <- rownames(geography.distance.df)
geography.distance.df.long <- melt(geography.distance.df, id.vars = "sampleID")
geography.distance.df.long <- na.omit(geography.distance.df.long)
climate.distance.mat <- as.matrix(climate.distance)
climate.distance.mat[upper.tri(climate.distance.mat, diag = TRUE)] <- NA
climate.distance.df <- as.data.frame(climate.distance.mat)
climate.distance.df$sampleID <- rownames(climate.distance.df)
climate.distance.df.long <- melt(climate.distance.df, id.vars = "sampleID")
climate.distance.df.long <- na.omit(climate.distance.df.long)
diet.distance.mat <- as.matrix(diet.distance)
diet.distance.mat[upper.tri(diet.distance.mat, diag = TRUE)] <- NA
diet.distance.df <- as.data.frame(diet.distance.mat)
diet.distance.df$sampleID <- rownames(diet.distance.df)
diet.distance.df.long <- melt(diet.distance.df, id.vars = "sampleID")
diet.distance.df.long <- na.omit(diet.distance.df.long)
names(bacteria.distance.df.long) <- c("sampleID","variable","BC")
bacteria.distance.df.long$subspecies <- "all"
bacteria.distance.df.long$Geography <- geography.distance.df.long$value
bacteria.distance.df.long$Climate <- climate.distance.df.long$value
bacteria.distance.df.long$Diet <- diet.distance.df.long$value
fig5all <- melt(bacteria.distance.df.long,
                id.vars = c("sampleID","variable","BC","subspecies"),
                measure.vars = c("Geography","Climate","Diet"))
names(fig5all) <- c("sampleID","variable","BC","subspecies","dataset","dist")

# Schwein
schwein.bacteria.distance.mat <- as.matrix(schwein.bacteria.distance)
schwein.bacteria.distance.mat[upper.tri(schwein.bacteria.distance.mat, diag = TRUE)] <- NA
schwein.bacteria.distance.df <- as.data.frame(schwein.bacteria.distance.mat)
schwein.bacteria.distance.df$sampleID <- rownames(schwein.bacteria.distance.df)
schwein.bacteria.distance.df.long <- melt(schwein.bacteria.distance.df, id.vars = "sampleID")
schwein.bacteria.distance.df.long <- na.omit(schwein.bacteria.distance.df.long)
schwein.geography.distance.mat[upper.tri(schwein.geography.distance.mat, diag = TRUE)] <- NA
schwein.geography.distance.df <- as.data.frame(schwein.geography.distance.mat)
schwein.geography.distance.df$sampleID <- rownames(schwein.geography.distance.df)
schwein.geography.distance.df.long <- melt(schwein.geography.distance.df, id.vars = "sampleID")
schwein.geography.distance.df.long <- na.omit(schwein.geography.distance.df.long)
schwein.climate.distance.mat <- as.matrix(schwein.climate.distance)
schwein.climate.distance.mat[upper.tri(schwein.climate.distance.mat, diag = TRUE)] <- NA
schwein.climate.distance.df <- as.data.frame(schwein.climate.distance.mat)
schwein.climate.distance.df$sampleID <- rownames(schwein.climate.distance.df)
schwein.climate.distance.df.long <- melt(schwein.climate.distance.df, id.vars = "sampleID")
schwein.climate.distance.df.long <- na.omit(schwein.climate.distance.df.long)
schwein.diet.distance.mat <- as.matrix(schwein.diet.distance)
schwein.diet.distance.mat[upper.tri(schwein.diet.distance.mat, diag = TRUE)] <- NA
schwein.diet.distance.df <- as.data.frame(schwein.diet.distance.mat)
schwein.diet.distance.df$sampleID <- rownames(schwein.diet.distance.df)
schwein.diet.distance.df.long <- melt(schwein.diet.distance.df, id.vars = "sampleID")
schwein.diet.distance.df.long <- na.omit(schwein.diet.distance.df.long)
names(schwein.bacteria.distance.df.long) <- c("sampleID","variable","BC")
schwein.bacteria.distance.df.long$subspecies <- "schweinfurthii"
schwein.bacteria.distance.df.long$Geography <- schwein.geography.distance.df.long$value
schwein.bacteria.distance.df.long$Climate <- schwein.climate.distance.df.long$value
schwein.bacteria.distance.df.long$Diet <- schwein.diet.distance.df.long$value
fig5schwein <- melt(schwein.bacteria.distance.df.long,
                id.vars = c("sampleID","variable","BC","subspecies"),
                measure.vars = c("Geography","Climate","Diet"))
names(fig5schwein) <- c("sampleID","variable","BC","subspecies","dataset","dist")

# Trog
trog.bacteria.distance.mat <- as.matrix(trog.bacteria.distance)
trog.bacteria.distance.mat[upper.tri(trog.bacteria.distance.mat, diag = TRUE)] <- NA
trog.bacteria.distance.df <- as.data.frame(trog.bacteria.distance.mat)
trog.bacteria.distance.df$sampleID <- rownames(trog.bacteria.distance.df)
trog.bacteria.distance.df.long <- melt(trog.bacteria.distance.df, id.vars = "sampleID")
trog.bacteria.distance.df.long <- na.omit(trog.bacteria.distance.df.long)
trog.geography.distance.mat[upper.tri(trog.geography.distance.mat, diag = TRUE)] <- NA
trog.geography.distance.df <- as.data.frame(trog.geography.distance.mat)
trog.geography.distance.df$sampleID <- rownames(trog.geography.distance.df)
trog.geography.distance.df.long <- melt(trog.geography.distance.df, id.vars = "sampleID")
trog.geography.distance.df.long <- na.omit(trog.geography.distance.df.long)
trog.climate.distance.mat <- as.matrix(trog.climate.distance)
trog.climate.distance.mat[upper.tri(trog.climate.distance.mat, diag = TRUE)] <- NA
trog.climate.distance.df <- as.data.frame(trog.climate.distance.mat)
trog.climate.distance.df$sampleID <- rownames(trog.climate.distance.df)
trog.climate.distance.df.long <- melt(trog.climate.distance.df, id.vars = "sampleID")
trog.climate.distance.df.long <- na.omit(trog.climate.distance.df.long)
trog.diet.distance.mat <- as.matrix(trog.diet.distance)
trog.diet.distance.mat[upper.tri(trog.diet.distance.mat, diag = TRUE)] <- NA
trog.diet.distance.df <- as.data.frame(trog.diet.distance.mat)
trog.diet.distance.df$sampleID <- rownames(trog.diet.distance.df)
trog.diet.distance.df.long <- melt(trog.diet.distance.df, id.vars = "sampleID")
trog.diet.distance.df.long <- na.omit(trog.diet.distance.df.long)
names(trog.bacteria.distance.df.long) <- c("sampleID","variable","BC")
trog.bacteria.distance.df.long$subspecies <- "troglodytes"
trog.bacteria.distance.df.long$Geography <- trog.geography.distance.df.long$value
trog.bacteria.distance.df.long$Climate <- trog.climate.distance.df.long$value
trog.bacteria.distance.df.long$Diet <- trog.diet.distance.df.long$value
fig5trog <- melt(trog.bacteria.distance.df.long,
                    id.vars = c("sampleID","variable","BC","subspecies"),
                    measure.vars = c("Geography","Climate","Diet"))
names(fig5trog) <- c("sampleID","variable","BC","subspecies","dataset","dist")

# Verus
verus.bacteria.distance.mat <- as.matrix(verus.bacteria.distance)
verus.bacteria.distance.mat[upper.tri(verus.bacteria.distance.mat, diag = TRUE)] <- NA
verus.bacteria.distance.df <- as.data.frame(verus.bacteria.distance.mat)
verus.bacteria.distance.df$sampleID <- rownames(verus.bacteria.distance.df)
verus.bacteria.distance.df.long <- melt(verus.bacteria.distance.df, id.vars = "sampleID")
verus.bacteria.distance.df.long <- na.omit(verus.bacteria.distance.df.long)
verus.geography.distance.mat[upper.tri(verus.geography.distance.mat, diag = TRUE)] <- NA
verus.geography.distance.df <- as.data.frame(verus.geography.distance.mat)
verus.geography.distance.df$sampleID <- rownames(verus.geography.distance.df)
verus.geography.distance.df.long <- melt(verus.geography.distance.df, id.vars = "sampleID")
verus.geography.distance.df.long <- na.omit(verus.geography.distance.df.long)
verus.climate.distance.mat <- as.matrix(verus.climate.distance)
verus.climate.distance.mat[upper.tri(verus.climate.distance.mat, diag = TRUE)] <- NA
verus.climate.distance.df <- as.data.frame(verus.climate.distance.mat)
verus.climate.distance.df$sampleID <- rownames(verus.climate.distance.df)
verus.climate.distance.df.long <- melt(verus.climate.distance.df, id.vars = "sampleID")
verus.climate.distance.df.long <- na.omit(verus.climate.distance.df.long)
verus.diet.distance.mat <- as.matrix(verus.diet.distance)
verus.diet.distance.mat[upper.tri(verus.diet.distance.mat, diag = TRUE)] <- NA
verus.diet.distance.df <- as.data.frame(verus.diet.distance.mat)
verus.diet.distance.df$sampleID <- rownames(verus.diet.distance.df)
verus.diet.distance.df.long <- melt(verus.diet.distance.df, id.vars = "sampleID")
verus.diet.distance.df.long <- na.omit(verus.diet.distance.df.long)
names(verus.bacteria.distance.df.long) <- c("sampleID","variable","BC")
verus.bacteria.distance.df.long$subspecies <- "verus"
verus.bacteria.distance.df.long$Geography <- verus.geography.distance.df.long$value
verus.bacteria.distance.df.long$Climate <- verus.climate.distance.df.long$value
verus.bacteria.distance.df.long$Diet <- verus.diet.distance.df.long$value
fig5verus <- melt(verus.bacteria.distance.df.long,
                 id.vars = c("sampleID","variable","BC","subspecies"),
                 measure.vars = c("Geography","Climate","Diet"))
names(fig5verus) <- c("sampleID","variable","BC","subspecies","dataset","dist")
fig5comb <- rbind(fig5all, fig5schwein, fig5trog, fig5verus)
fig5comb$subspecies <- as.factor(fig5comb$subspecies)

# Figure
# Figure out midpoint
(max(bacteria.distance.df.long$Geography) - min(bacteria.distance.df.long$Geography))/2
(max(bacteria.distance.df.long$Climate) - min(bacteria.distance.df.long$Climate))/2
(max(bacteria.distance.df.long$Diet) - min(bacteria.distance.df.long$Diet))/2
label.df.5 <- data.frame(subspecies = c("all","all","all",
                                        "schweinfurthii","schweinfurthii","schweinfurthii",
                                        "troglodytes","troglodytes","troglodytes",
                                        "verus","verus","verus"),
                         dataset = c("Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet"),
                         x = c(23.86629, 2.93386, 0.5,
                               23.86629, 2.93386, 0.5,
                               23.86629, 2.93386, 0.5,
                               23.86629, 2.93386, 0.5),
                         y = c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25, 0.25, 0.25),
                         label = c("n = 560, r = 0.68, p = 0.001",
                                   "n = 560, r = 0.40, p = 0.001",
                                   "n = 560, r = 0.33, p = 0.001",
                                   "n = 134, r = 0.51, p = 0.001",
                                   "n = 134, r = 0.58, p = 0.001",
                                   "n = 134, r = 0.60, p = 0.001",
                                   "n = 86, r = 0.52, p = 0.001",
                                   "n = 86, r = 0.53, p = 0.001",
                                   "n = 86, r = 0.30, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001",
                                   "n = 312, r = 0.24, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001"))
ggplot(data = fig5comb, aes(dist, BC)) +
  geom_point(size = 0.5, alpha = 0.01) +
  geom_smooth(data = subset(fig5comb, dataset == "Geography"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5comb, dataset == "Climate"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5comb, dataset == "Diet"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_text(data = label.df.5, aes(x = x, y = y, label = label), size = 2) +
  labs(x = "Distance",
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x") +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))



##################### __Part 2: Add humans ####################################################
humans <- select(humans, -X)
names(humans) == names(fig5comb)
fig5whuman <- rbind(fig5comb, humans)
(max(subset(fig5whuman, dataset == "Geography")$dist) - min(subset(fig5whuman, dataset == "Geography")$dist))/2
(max(subset(fig5whuman, dataset == "Climate")$dist) - min(subset(fig5whuman, dataset == "Climate")$dist))/2

label.df.h <- data.frame(subspecies = c("all","all","all",
                                        "schweinfurthii","schweinfurthii","schweinfurthii",
                                        "troglodytes","troglodytes","troglodytes",
                                        "verus","verus","verus",
                                        "human","human"),
                         dataset = c("Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate","Diet",
                                     "Geography","Climate"),
                         x = c(28.1177, 3.602599, 0.5,
                               28.1177, 3.602599, 0.5,
                               28.1177, 3.602599, 0.5,
                               28.1177, 3.602599, 0.5,
                               28.1177, 3.602599),
                         y = c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125, 0.125,
                               0.125,0.125, 0.125),
                         label = c("n = 560, r = 0.68, p = 0.001",
                                   "n = 560, r = 0.40, p = 0.001",
                                   "n = 560, r = 0.33, p = 0.001",
                                   "n = 134, r = 0.51, p = 0.001",
                                   "n = 134, r = 0.58, p = 0.001",
                                   "n = 134, r = 0.60, p = 0.001",
                                   "n = 86, r = 0.52, p = 0.001",
                                   "n = 86, r = 0.53, p = 0.001",
                                   "n = 86, r = 0.30, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001",
                                   "n = 312, r = 0.24, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001",
                                   "n = 2296, r = 0.007, p = 0.05",
                                   "n = 2296, r = 0.007, p = 0.29"))
facet_names_5 <- c("Geography" = "Geography", "Climate" = "Climate", "Diet" = "Diet/Tool Use",
                   "all" = "all", "schweinfurthii" = "schweinfurthii", 
                   "troglodytes" = "troglodytes", "verus" = "verus", "human" = "human")
ggplot(data = fig5whuman, aes(dist, BC)) +
  geom_point(data = subset(fig5whuman, subspecies != "human"),
             size = 0.5, alpha = 0.01) +
  geom_point(data = subset(fig5whuman, subspecies == "human"),
             size = 0.1, alpha = 0.01) +
  geom_smooth(data = subset(fig5whuman, dataset == "Geography" & subspecies != "human"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5whuman, dataset == "Diet"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whuman, dataset == "Climate" & subspecies == "all"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whuman, dataset == "Climate" & subspecies == "schweinfurthii"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whuman, dataset == "Climate" & subspecies == "troglodytes"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whuman, dataset == "Climate" & subspecies == "verus"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whuman, subspecies == "human"),
              method = lm, formula = y ~ x, size = 0.5, linetype = "dashed") +
  geom_text(data = label.df.h, aes(x = x, y = y, label = label), size = 2) +
  labs(x = "Distance",
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x", labeller = as_labeller(facet_names_5)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))



##################### __Part 3: Add Parasites #################################################
# Could also add another column with parasite Jaccard!
# Take object bcj from section 5S
# Update: We've decided to keep this as a separate figure
bcjtocomb <- bcj
names(bcjtocomb) <- c("sampleID","variable","BC","dist","subspecies")
bcjtocomb$dataset <- "Parasites"
bcjtocomb <- select(bcjtocomb, sampleID, variable, BC, subspecies, dataset, dist)
names(bcjtocomb) == names(fig5whuman)
fig5whumanwpara <- rbind(fig5whuman, bcjtocomb)
fig5whumanwpara$dataset <- factor(fig5whumanwpara$dataset,
                                  levels = c("Parasites","Diet","Climate","Geography"))
levels(fig5whumanwpara$dataset)
fig5whumanwpara <- subset(fig5whumanwpara, subspecies != "ellioti")
label.df.p <- data.frame(subspecies = c("all","all","all","all",
                                        "schweinfurthii","schweinfurthii","schweinfurthii","schweinfurthii",
                                        "troglodytes","troglodytes","troglodytes","troglodytes",
                                        "verus","verus","verus","verus",
                                        "human","human"),
                         dataset = c("Parasites","Geography","Climate","Diet",
                                     "Parasites","Geography","Climate","Diet",
                                     "Parasites","Geography","Climate","Diet",
                                     "Parasites","Geography","Climate","Diet",
                                     "Geography","Climate"),
                         x = c(0.5, 28.1177, 3.602599, 0.5,
                               0.5, 28.1177, 3.602599, 0.5,
                               0.5, 28.1177, 3.602599, 0.5,
                               0.5, 28.1177, 3.602599, 0.5,
                               28.1177, 3.602599),
                         y = c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,
                               0.125,0.125,0.125,0.125, 0.125,0.125,0.125, 0.125),
                         label = c("n = 560, r = 0.59, p = 0.001",
                                   "n = 560, r = 0.68, p = 0.001",
                                   "n = 560, r = 0.40, p = 0.001",
                                   "n = 560, r = 0.33, p = 0.001",
                                   "n = 134, r = 0.36, p = 0.001",
                                   "n = 134, r = 0.51, p = 0.001",
                                   "n = 134, r = 0.58, p = 0.001",
                                   "n = 134, r = 0.60, p = 0.001",
                                   "n = 86, r = 0.48, p = 0.001",
                                   "n = 86, r = 0.52, p = 0.001",
                                   "n = 86, r = 0.53, p = 0.001",
                                   "n = 86, r = 0.30, p = 0.001",
                                   "n = 312, r = 0.34, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001",
                                   "n = 312, r = 0.24, p = 0.001",
                                   "n = 312, r = 0.29, p = 0.001",
                                   "n = 2296, r = 0.007, p = 0.05",
                                   "n = 2296, r = 0.007, p = 0.29"))

png(file = "forPPTFigure5whumanwparaHex.png", width = 7, height = 7, unit = "in", res = 300)
ggplot(data = fig5whumanwpara, aes(dist, BC)) +
  geom_point(data = subset(fig5whumanwpara, subspecies != "human"), size = 0.5, alpha = 0.01) +
  geom_hex(data = subset(fig5whumanwpara, subspecies == "human"), aes(fill = stat(count)),
           bins = 18) +
  scale_fill_gradientn(colors = get.palette(11),na.value = 'white',trans = 'sqrt', 
                       breaks = pretty_breaks(n = 4)) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Diet" | dataset == "Parasites"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Climate" & subspecies == "all"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Climate" & subspecies == "schweinfurthii"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Climate" & subspecies == "troglodytes"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Climate" & subspecies == "verus"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Geography" & subspecies == "all"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Geography" & subspecies == "schweinfurthii"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Geography" & subspecies == "troglodytes"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5whumanwpara, dataset == "Geography" & subspecies == "verus"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_text(data = label.df.p, aes(x = x, y = y, label = label), size = 2) +
  labs(x = "Distance",
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x") +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))
dev.off()



#### __Update: Add vegetation and remove climate and diet ####
mantel(veg.distance, bacteria.distance.v) # r = 0.54, p = 0.001
mantel(elli.veg.distance, elli.bacteria.distance.v) # r = 0.66, p = 0.001
mantel(schwein.veg.distance, schwein.bacteria.distance.v) # r = 0.61, p = 0.001
mantel(trog.veg.distance, trog.bacteria.distance.v) # r = 0.52, p = 0.001
mantel(verus.veg.distance, verus.bacteria.distance.v) # r = 0.29, p = 0.001

mantel.partial(veg.distance, bacteria.distance.v, geog.distance.v) # r = 0.27, p = 0.001
mantel.partial(elli.veg.distance,elli.bacteria.distance.v,elli.geog.distance.v) # NA
mantel.partial(schwein.veg.distance,schwein.bacteria.distance.v,schwein.geog.distance.v)#.4,.001
mantel.partial(trog.veg.distance, trog.bacteria.distance.v, trog.geog.distance.v) # 0.29, 0.001
mantel.partial(verus.veg.distance, verus.bacteria.distance.v, verus.geog.distance.v) #0.09,0.001

# Assemble Figure Dataframe. Geography and Diet are the same. Use subset bac, veg instead of clim
bacteria.distance.v.mat <- as.matrix(bacteria.distance.v)
bacteria.distance.v.mat[upper.tri(bacteria.distance.v.mat, diag = TRUE)] <- NA
bacteria.distance.v.df <- as.data.frame(bacteria.distance.v.mat)
bacteria.distance.v.df$sampleID <- rownames(bacteria.distance.v.df)
bacteria.distance.v.df.long <- melt(bacteria.distance.v.df, id.vars = "sampleID")
bacteria.distance.v.df.long <- na.omit(bacteria.distance.v.df.long)
veg.distance.mat <- as.matrix(veg.distance)
veg.distance.mat[upper.tri(veg.distance.mat, diag = TRUE)] <- NA
veg.distance.df <- as.data.frame(veg.distance.mat)
veg.distance.df$sampleID <- rownames(veg.distance.df)
veg.distance.df.long <- melt(veg.distance.df, id.vars = "sampleID")
veg.distance.df.long <- na.omit(veg.distance.df.long)
names(bacteria.distance.v.df.long) <- c("sampleID","variable","BC")
bacteria.distance.v.df.long$subspecies <- "all"
bacteria.distance.v.df.long$Veg <- veg.distance.df.long$value
fig5all.vGD <- subset(fig5all, dataset != "Climate")
fig5all.vV <- melt(bacteria.distance.v.df.long,
                   id.vars = c("sampleID","variable","BC","subspecies"),
                   measure.vars = c("Veg"))
names(fig5all.vV) <- c("sampleID","variable","BC","subspecies","dataset","dist")
fig5all.v <- rbind(fig5all.vGD, fig5all.vV)

# Schwein
schwein.bacteria.distance.v.mat <- as.matrix(schwein.bacteria.distance.v)
schwein.bacteria.distance.v.mat[upper.tri(schwein.bacteria.distance.v.mat, diag = TRUE)] <- NA
schwein.bacteria.distance.v.df <- as.data.frame(schwein.bacteria.distance.v.mat)
schwein.bacteria.distance.v.df$sampleID <- rownames(schwein.bacteria.distance.v.df)
schwein.bacteria.distance.v.df.long <- melt(schwein.bacteria.distance.v.df, id.vars = "sampleID")
schwein.bacteria.distance.v.df.long <- na.omit(schwein.bacteria.distance.v.df.long)
schwein.veg.distance.mat <- as.matrix(schwein.veg.distance)
schwein.veg.distance.mat[upper.tri(schwein.veg.distance.mat, diag = TRUE)] <- NA
schwein.veg.distance.df <- as.data.frame(schwein.veg.distance.mat)
schwein.veg.distance.df$sampleID <- rownames(schwein.veg.distance.df)
schwein.veg.distance.df.long <- melt(schwein.veg.distance.df, id.vars = "sampleID")
schwein.veg.distance.df.long <- na.omit(schwein.veg.distance.df.long)
names(schwein.bacteria.distance.v.df.long) <- c("sampleID","variable","BC")
schwein.bacteria.distance.v.df.long$subspecies <- "schweinfurthii"
schwein.bacteria.distance.v.df.long$Veg <- schwein.veg.distance.df.long$value
schwein.fig5all.vGD <- subset(fig5schwein, dataset != "Climate")
schwein.fig5all.vV <- melt(schwein.bacteria.distance.v.df.long,
                           id.vars = c("sampleID","variable","BC","subspecies"),
                           measure.vars = c("Veg"))
names(schwein.fig5all.vV) <- c("sampleID","variable","BC","subspecies","dataset","dist")
schwein.fig5all.v <- rbind(schwein.fig5all.vGD, schwein.fig5all.vV)

# Trog
trog.bacteria.distance.v.mat <- as.matrix(trog.bacteria.distance.v)
trog.bacteria.distance.v.mat[upper.tri(trog.bacteria.distance.v.mat, diag = TRUE)] <- NA
trog.bacteria.distance.v.df <- as.data.frame(trog.bacteria.distance.v.mat)
trog.bacteria.distance.v.df$sampleID <- rownames(trog.bacteria.distance.v.df)
trog.bacteria.distance.v.df.long <- melt(trog.bacteria.distance.v.df, id.vars = "sampleID")
trog.bacteria.distance.v.df.long <- na.omit(trog.bacteria.distance.v.df.long)
trog.veg.distance.mat <- as.matrix(trog.veg.distance)
trog.veg.distance.mat[upper.tri(trog.veg.distance.mat, diag = TRUE)] <- NA
trog.veg.distance.df <- as.data.frame(trog.veg.distance.mat)
trog.veg.distance.df$sampleID <- rownames(trog.veg.distance.df)
trog.veg.distance.df.long <- melt(trog.veg.distance.df, id.vars = "sampleID")
trog.veg.distance.df.long <- na.omit(trog.veg.distance.df.long)
names(trog.bacteria.distance.v.df.long) <- c("sampleID","variable","BC")
trog.bacteria.distance.v.df.long$subspecies <- "troglodytes"
trog.bacteria.distance.v.df.long$Veg <- trog.veg.distance.df.long$value
trog.fig5all.vGD <- subset(fig5trog, dataset != "Climate")
trog.fig5all.vV <- melt(trog.bacteria.distance.v.df.long,
                        id.vars = c("sampleID","variable","BC","subspecies"),
                        measure.vars = c("Veg"))
names(trog.fig5all.vV) <- c("sampleID","variable","BC","subspecies","dataset","dist")
trog.fig5all.v <- rbind(trog.fig5all.vGD, trog.fig5all.vV)

# Verus
verus.bacteria.distance.v.mat <- as.matrix(verus.bacteria.distance.v)
verus.bacteria.distance.v.mat[upper.tri(verus.bacteria.distance.v.mat, diag = TRUE)] <- NA
verus.bacteria.distance.v.df <- as.data.frame(verus.bacteria.distance.v.mat)
verus.bacteria.distance.v.df$sampleID <- rownames(verus.bacteria.distance.v.df)
verus.bacteria.distance.v.df.long <- melt(verus.bacteria.distance.v.df, id.vars = "sampleID")
verus.bacteria.distance.v.df.long <- na.omit(verus.bacteria.distance.v.df.long)
verus.veg.distance.mat <- as.matrix(verus.veg.distance)
verus.veg.distance.mat[upper.tri(verus.veg.distance.mat, diag = TRUE)] <- NA
verus.veg.distance.df <- as.data.frame(verus.veg.distance.mat)
verus.veg.distance.df$sampleID <- rownames(verus.veg.distance.df)
verus.veg.distance.df.long <- melt(verus.veg.distance.df, id.vars = "sampleID")
verus.veg.distance.df.long <- na.omit(verus.veg.distance.df.long)
names(verus.bacteria.distance.v.df.long) <- c("sampleID","variable","BC")
verus.bacteria.distance.v.df.long$subspecies <- "verus"
verus.bacteria.distance.v.df.long$Veg <- verus.veg.distance.df.long$value
verus.fig5all.vGD <- subset(fig5verus, dataset != "Climate")
verus.fig5all.vV <- melt(verus.bacteria.distance.v.df.long,
                         id.vars = c("sampleID","variable","BC","subspecies"),
                         measure.vars = c("Veg"))
names(verus.fig5all.vV) <- c("sampleID","variable","BC","subspecies","dataset","dist")
verus.fig5all.v <- rbind(verus.fig5all.vGD, verus.fig5all.vV)
fig5comb.v <- rbind(fig5all.v, schwein.fig5all.v, trog.fig5all.v, verus.fig5all.v)
fig5comb.v$subspecies <- as.factor(fig5comb.v$subspecies)
levels(fig5comb.v$subspecies)
levels(fig5comb.v$dataset)
fig5comb.v$dataset <- factor(fig5comb.v$dataset, levels = c("Geography","Veg","Diet"))

# Figure
(max(bacteria.distance.v.df.long$Veg) - min(bacteria.distance.v.df.long$Veg))/2
label.df.5.v <- data.frame(subspecies = c("all","all","all",
                                          "schweinfurthii","schweinfurthii","schweinfurthii",
                                          "troglodytes","troglodytes","troglodytes",
                                          "verus","verus","verus"),
                           dataset = c("Geography","Veg","Diet",
                                       "Geography","Veg","Diet",
                                       "Geography","Veg","Diet",
                                       "Geography","Veg","Diet"),
                           x = c(23.86629, 0.5, 0.5,
                                 23.86629, 0.5, 0.5,
                                 23.86629, 0.5, 0.5,
                                 23.86629, 0.5, 0.5),
                           y = c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25, 0.25, 0.25),
                           label = c("n = 560, r = 0.68, p = 0.001",
                                     "n = 509, r = 0.54, p = 0.001",
                                     "n = 560, r = 0.33, p = 0.001",
                                     "n = 134, r = 0.51, p = 0.001",
                                     "n = 134, r = 0.61, p = 0.001",
                                     "n = 134, r = 0.60, p = 0.001",
                                     "n = 86, r = 0.52, p = 0.001",
                                     "n = 35, r = 0.52, p = 0.001",
                                     "n = 86, r = 0.30, p = 0.001",
                                     "n = 312, r = 0.29, p = 0.001",
                                     "n = 312, r = 0.29, p = 0.001",
                                     "n = 312, r = 0.29, p = 0.001"))
facet_names_5.v <- c("Geography" = "Geography", "Veg" = "Vegetation", "Diet" = "Tool Use",
                     "all" = "all", "schweinfurthii" = "schweinfurthii", 
                     "troglodytes" = "troglodytes", "verus" = "verus", "human" = "human")

# Update to remove diet
fig5comb.v.nd <- subset(fig5comb.v, dataset != "Diet")
label.df.5.v <- data.frame(subspecies = c("all","all",
                                          "schweinfurthii","schweinfurthii",
                                          "troglodytes","troglodytes",
                                          "verus","verus"),
                           dataset = c("Geography","Veg",
                                       "Geography","Veg",
                                       "Geography","Veg",
                                       "Geography","Veg"),
                           x = c(2650, 0.5,
                                 2650, 0.5,
                                 2650, 0.5,
                                 2650, 0.5),
                           y = c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25),
                           label = c("n = 560, r = 0.68, p = 0.001",
                                     "n = 509, r = 0.54, p = 0.001",
                                     "n = 134, r = 0.51, p = 0.001",
                                     "n = 134, r = 0.61, p = 0.001",
                                     "n = 86, r = 0.52, p = 0.001",
                                     "n = 35, r = 0.52, p = 0.001",
                                     "n = 312, r = 0.29, p = 0.001",
                                     "n = 312, r = 0.29, p = 0.001"))
facet_names_5.v <- c("Geography" = "a) Geography", "Veg" = "b) Vegetation",
                     "all" = "all", "schweinfurthii" = "East", 
                     "troglodytes" = "Central", "verus" = "West")
levels(fig5comb.v.nd$subspecies)
fig5comb.v.nd$subspecies <- factor(fig5comb.v.nd$subspecies,
                                   levels = c("all", "verus", "troglodytes", "schweinfurthii"))
# pdf(file = "Figure5.pdf", width = 5, height = 5)
tiff(file = "PNGs/Figure4.tiff", width = 5, height = 5, units = "in", res = 300)
ggplot(data = fig5comb.v.nd, aes(dist, BC)) +
  geom_point(size = 0.5, alpha = 0.01) +
  geom_smooth(data = subset(fig5comb.v.nd, dataset == "Geography"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig5comb.v.nd, dataset == "Veg"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_text(data = label.df.5.v, aes(x = x, y = y, label = label), size = 3) +
  labs(x = "Distance",
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x", labeller = as_labeller(facet_names_5.v)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10),
        panel.spacing.x = unit(1, "lines"))
dev.off()




##################### _6. Multipanel - BC and Jaccard within vs. between sites ################
##################### __Part 1: Bray Curtis ###################################################
# 5 panels, all and each subsp.
# Make long dataframe, add Site and subspecies
# Within versus among sites, all samples
# Within versus among sites within each subspecies
bac_bray_mat <- as.matrix(bacteria.distance)
bac_bray_mat[upper.tri(bac_bray_mat, diag = TRUE)] <- NA
bac_bray_df <- as.data.frame(bac_bray_mat)
bac_bray_df$sampleID <- rownames(bac_bray_df)
bac_bray_df_long <- melt(bac_bray_df, id.vars = "sampleID")
bac_bray_df_long <- na.omit(bac_bray_df_long)
bac_bray_df_long$sampleID <- as.factor(bac_bray_df_long$sampleID)
# Now add site and subspecies, matching to col1 and col2
site_subspecies <- select(bac_euk$map_loaded, sampleID, Site, subspecies)
bac_bray_df_long <- inner_join(bac_bray_df_long, site_subspecies, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long <- inner_join(bac_bray_df_long, site_subspecies, 
                               by = c("variable" = "sampleID"))
# Make new column indicating if comparison is within or between sites
for (i in 1:nrow(bac_bray_df_long)) {
  ifelse(bac_bray_df_long$Site.x[i] == bac_bray_df_long$Site.y[i],
         bac_bray_df_long$comparison[i] <- "within",
         bac_bray_df_long$comparison[i] <- "between")
}
bac_bray_df_long$comparison <- as.factor(bac_bray_df_long$comparison)
# Subspecies dfs
bac_bray_df_long_elli <- subset(bac_bray_df_long, 
                                subspecies.x == "ellioti" & subspecies.y == "ellioti")
bac_bray_df_long_schwein <- subset(bac_bray_df_long, 
                                   subspecies.x=="schweinfurthii"&subspecies.y=="schweinfurthii")
bac_bray_df_long_trog <- subset(bac_bray_df_long, 
                                subspecies.x == "troglodytes" & subspecies.y == "troglodytes")
bac_bray_df_long_verus <- subset(bac_bray_df_long, subspecies.x == "verus" & subspecies.y == "verus")
# Sample sizes
table(bac_bray_df_long$comparison)
table(bac_bray_df_long_elli$comparison)
table(bac_bray_df_long_schwein$comparison)
table(bac_bray_df_long_trog$comparison)
table(bac_bray_df_long_verus$comparison)
# Stats (all highly significant)
t.test(value ~ comparison, data = bac_bray_df_long)
t.test(value ~ comparison, data = bac_bray_df_long_elli)
t.test(value ~ comparison, data = bac_bray_df_long_schwein)
t.test(value ~ comparison, data = bac_bray_df_long_trog)
t.test(value ~ comparison, data = bac_bray_df_long_verus)
wilcox.test(value ~ comparison, data = bac_bray_df_long)
wilcox.test(value ~ comparison, data = bac_bray_df_long_elli)
wilcox.test(value ~ comparison, data = bac_bray_df_long_schwein)
wilcox.test(value ~ comparison, data = bac_bray_df_long_trog)
wilcox.test(value ~ comparison, data = bac_bray_df_long_verus)
# Combine for multipanel
bac_bray_df_long$subspecies <- "all"
bac_bray_df_long_elli$subspecies <- "ellioti"
bac_bray_df_long_schwein$subspecies <- "schweinfurthii"
bac_bray_df_long_trog$subspecies <- "troglodytes"
bac_bray_df_long_verus$subspecies <- "verus"
bac_bray_df_long_multi <- rbind(bac_bray_df_long,bac_bray_df_long_elli,bac_bray_df_long_schwein,
                                bac_bray_df_long_trog,bac_bray_df_long_verus)
bac_bray_df_long_multi$subspecies <- as.factor(bac_bray_df_long_multi$subspecies)
levels(bac_bray_df_long_multi$subspecies)
bac_bray_df_long_multi$subspecies <- factor(bac_bray_df_long_multi$subspecies,
                                levels = c("all","verus","ellioti","troglodytes","schweinfurthii"))
# Graph
facet.names.6a <- c('all' = "a) all",
                    'ellioti' = "c) N-C",
                    'schweinfurthii' = "e) East",
                    'troglodytes' = "d) Central",
                    'verus' = "b) West")
label.df.6a <- data.frame(subspecies = c("all","all","ellioti","ellioti",
                                         "schweinfurthii","schweinfurthii",
                                         "troglodytes","troglodytes","verus","verus"),
                          comparison = c("between","within","between","within","between","within",
                                         "between","within","between","within"),
                          Value = c(1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05),
                          Sig = c("a","b","a","b","a","b","a","b","a","b"))
ggplot(data = bac_bray_df_long_multi, aes(comparison, value)) +
  geom_jitter(data = subset(bac_bray_df_long_multi, subspecies == "all" | subspecies == "verus"),
              size = 0.5, alpha = 0.01) +
  geom_jitter(data = subset(bac_bray_df_long_multi, subspecies == "ellioti"),
              size = 0.5, alpha = 0.06) +
  geom_jitter(data = subset(bac_bray_df_long_multi, subspecies == "schweinfurthii" |
                              subspecies == "troglodytes"), size = 0.5, alpha = 0.02) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "blue") +
  geom_text(data = label.df.6a, aes(x = comparison, y = Value, label = Sig, group=NULL)) +
  labs(x = "Site Comparison",
       y = "Bray-Curtis Dissimilarity") +
  facet_wrap(~ subspecies, ncol = 5, labeller = as_labeller(facet.names.6a)) +
  ylim(0,1.05) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 11))



####################### __Part 2: Jaccard #####################################################
# 5 panels, all and each subsp.
# Make long dataframe, add Site and subspecies
# Within versus among sites, all samples
# Within versus among sites within each subspecies
euk_jac_mat <- as.matrix(parasites.distance)
euk_jac_mat[upper.tri(euk_jac_mat, diag = TRUE)] <- NA
euk_jac_df <- as.data.frame(euk_jac_mat)
euk_jac_df$sampleID <- rownames(euk_jac_df)
euk_jac_df_long <- melt(euk_jac_df, id.vars = "sampleID")
euk_jac_df_long <- na.omit(euk_jac_df_long)
euk_jac_df_long$sampleID <- as.factor(euk_jac_df_long$sampleID)
# Now add site and subspecies, matching to col1 and col2
site_subspecies <- select(bac_euk$map_loaded, sampleID, Site, subspecies)
euk_jac_df_long <- inner_join(euk_jac_df_long, site_subspecies, 
                              by = c("sampleID" = "sampleID"))
euk_jac_df_long <- inner_join(euk_jac_df_long, site_subspecies, 
                              by = c("variable" = "sampleID"))
# Make new column indicating if comparison is within or between sites
for (i in 1:nrow(euk_jac_df_long)) {
  ifelse(euk_jac_df_long$Site.x[i] == euk_jac_df_long$Site.y[i],
         euk_jac_df_long$comparison[i] <- "within",
         euk_jac_df_long$comparison[i] <- "between")
}
euk_jac_df_long$comparison <- as.factor(euk_jac_df_long$comparison)
# Subspecies dfs
euk_jac_df_long_elli <- subset(euk_jac_df_long, 
                               subspecies.x == "ellioti" & subspecies.y == "ellioti")
euk_jac_df_long_schwein <- subset(euk_jac_df_long, 
                                  subspecies.x=="schweinfurthii"&subspecies.y=="schweinfurthii")
euk_jac_df_long_trog <- subset(euk_jac_df_long, 
                               subspecies.x == "troglodytes" & subspecies.y == "troglodytes")
euk_jac_df_long_verus <- subset(euk_jac_df_long, subspecies.x == "verus" & subspecies.y == "verus")
# Sample sizes
table(euk_jac_df_long$comparison)
table(euk_jac_df_long_elli$comparison)
table(euk_jac_df_long_schwein$comparison)
table(euk_jac_df_long_trog$comparison)
table(euk_jac_df_long_verus$comparison)
# Stats (all highly significant)
t.test(value ~ comparison, data = euk_jac_df_long)
t.test(value ~ comparison, data = euk_jac_df_long_elli)
t.test(value ~ comparison, data = euk_jac_df_long_schwein)
t.test(value ~ comparison, data = euk_jac_df_long_trog)
t.test(value ~ comparison, data = euk_jac_df_long_verus)
# Combine for multipanel
euk_jac_df_long$subspecies <- "all"
euk_jac_df_long_elli$subspecies <- "ellioti"
euk_jac_df_long_schwein$subspecies <- "schweinfurthii"
euk_jac_df_long_trog$subspecies <- "troglodytes"
euk_jac_df_long_verus$subspecies <- "verus"
euk_jac_df_long_multi <- rbind(euk_jac_df_long,euk_jac_df_long_elli,euk_jac_df_long_schwein,
                               euk_jac_df_long_trog,euk_jac_df_long_verus)
euk_jac_df_long_multi$subspecies <- as.factor(bac_bray_df_long_multi$subspecies)
levels(euk_jac_df_long_multi$subspecies)
euk_jac_df_long_multi$subspecies <- factor(euk_jac_df_long_multi$subspecies,
                                levels = c("all","verus","ellioti","troglodytes","schweinfurthii"))
# Graph
facet.names.6b <- c('all' = "a) all",
                    'ellioti' = "c) N-C",
                    'schweinfurthii' = "e) East",
                    'troglodytes' = "d) Central",
                    'verus' = "b) West")
label.df.6b <- data.frame(subspecies = c("all","all","ellioti","ellioti",
                                         "schweinfurthii","schweinfurthii",
                                         "troglodytes","troglodytes","verus","verus"),
                          comparison = c("between","within","between","within","between","within",
                                         "between","within","between","within"),
                          Value = c(1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05,1.05),
                          Sig = c("a","b","a","b","a","b","a","b","a","b"))
ggplot(data = euk_jac_df_long_multi, aes(comparison, value)) +
  geom_jitter(data = subset(euk_jac_df_long_multi, subspecies == "all" | subspecies == "verus"),
              size = 0.5, alpha = 0.01) +
  geom_jitter(data = subset(euk_jac_df_long_multi, subspecies == "ellioti"),
              size = 0.5, alpha = 0.06) +
  geom_jitter(data = subset(euk_jac_df_long_multi, subspecies == "schweinfurthii" |
                              subspecies == "troglodytes"), size = 0.5, alpha = 0.02) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "blue") +
  geom_text(data = label.df.6b, aes(x = comparison, y = Value, label = Sig, group=NULL)) +
  labs(x = "Site Comparison",
       y = "Jaccard Dissimilarity") +
  facet_wrap(~ subspecies, ncol = 5, labeller = as_labeller(facet.names.6b)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 11))



##################### __Part 3: Combine #######################################################
# Run after completing jaccard section
wabac <- bac_bray_df_long_multi
wabac$dissim <- "Prokaryotes"
wabac$dissim <- as.factor(wabac$dissim)
wapar <- euk_jac_df_long_multi
wapar$dissim <- "Parasites"
wapar$dissim <- as.factor(wapar$dissim)
names(wabac) == names(wapar)
withamon_comb <- rbind(wabac, wapar)
withamon_comb <- subset(withamon_comb, subspecies != "all")
levels(withamon_comb$subspecies)
withamon_comb$subspecies <- factor(withamon_comb$subspecies,
                                levels = c("all","verus","ellioti","troglodytes","schweinfurthii"))
# Graph
facet.names.6 <- c('all' = "all",
                   'ellioti' = "N-C",
                   'schweinfurthii' = "East",
                   'troglodytes' = "Central",
                   'verus' = "West",
                   'Prokaryotes' = "Prokaryotes",
                   'Parasites' = "Parasites")
pdf(file = "Figure6.pdf", width = 7, height = 4)
ggplot(data = withamon_comb, aes(comparison, value)) +
  geom_jitter(alpha = 0.09, size = 0.5) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "blue") +
  labs(x = "Site Comparison",
       y = "Dissimilarity") +
  facet_grid(dissim ~ subspecies, labeller = as_labeller(facet.names.6)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 9))
dev.off()



################################### _7. Diet ##################################################
# New diet analysis - just algae, honey, nuts, termites
# Test all data and just direct
# Run full Permanova's
# To avoid errors, clean packages and reload a few
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(plyr)
library(tidyverse)
library(mctoolsr)
library(vegan)
library(indicspecies)
library(phyloseq)
library(ggtree)
library(BiodiversityR)
diet_update <- readRDS("bac_euk.rds")
diet_update$map_loaded <- left_join(diet_update$map_loaded, plants, by = c("Site" = "sites"))
rownames(diet_update$map_loaded) <- rownames(bac_euk$map_loaded)
f <- c("algae","honey","nuts","termites","algaeD","honeyD","nutsD","termitesD","Diet","DietD")
diet_update <- diet_update$map_loaded %>% 
  select(-ants, -fruit, -marrow, -meat, -'palm heart', -tubers, -water, -Diet) %>%
  mutate(algae = replace(algae,
                         Site == "Issa",
                         1)) %>%
  mutate(honey = replace(honey,
                         Site == "Tai_Ecotourism",
                         1)) %>%
  mutate(honey = replace(honey,
                         Site == "Tai_Recherche",
                         1)) %>%
  mutate(termites = replace(termites,
                            Site == "Tai_Ecotourism",
                            1)) %>%
  mutate(termites = replace(termites,
                            Site == "Tai_Recherche",
                            1)) %>%
  mutate(algaeD = algae) %>%
  mutate(algaeD = replace(algaeD, 
                          Site == "Bakoun" | Site == "Sangaredi" | Site == "Sobeya",
                          0)) %>%
  mutate(honeyD = honey) %>%
  mutate(honeyD = replace(honeyD, 
                          Site == "Bafing" | Site == "Boe" | Site == "Bwindi" | Site == "Gishwati",
                          0)) %>%
  mutate(nutsD = nuts) %>%
  mutate(nutsD = replace(nutsD, 
                         Site == "East_Nimba" | Site == "Sapo",
                         0)) %>%
  mutate(termitesD = termites) %>%
  mutate(termitesD = replace(termitesD, 
                             Site == "Kayan" | Site == "Mt_Cameroon" | Site == "Sobeya",
                             0)) %>%
  mutate(Diet = paste(algae, honey, nuts, termites, sep = "")) %>%
  mutate(DietD = paste(algaeD, honeyD, nutsD, termitesD, sep = "")) %>%
  mutate_at(f, factor) %>%
  mutate(Diet = recode(Diet,
                       "0000" = "None",
                       "0001" = "Termites",
                       "0010" = "Nuts",
                       "0100" = "Honey",
                       "0101" = "Honey, Termites",
                       "0111" = "Nuts, Honey, Termites",
                       "1000" = "Algae",
                       "1001" = "Algae, Termites",
                       "1101" = "Algae, Honey, Termites")) %>%
  mutate(DietD = recode(DietD,
                        "0000" = "None",
                        "0001" = "Termites",
                        "0010" = "Nuts",
                        "0100" = "Honey",
                        "0101" = "Honey, Termites",
                        "0111" = "Nuts, Honey, Termites",
                        "1000" = "Algae",
                        "1001" = "Algae, Termites"))
rownames(diet_update) <- rownames(bac_euk$map_loaded)
levels(diet_update$Diet)

#### __Part 1: Bray-Curtis ####
# a) ellioti
pcoa.bc.elli <- cmdscale(elli.bacteria.distance, k = nrow(elli$map_loaded)-1, eig=T)
eigenvals(pcoa.bc.elli)/sum(eigenvals(pcoa.bc.elli)) # 26.6, 10.3
elli$map_loaded$Axis01 <- scores(pcoa.bc.elli)[,1]
elli$map_loaded$Axis02 <- scores(pcoa.bc.elli)[,2]
# b) schweinfurthii
pcoa.bc.schwein <- cmdscale(schwein.bacteria.distance, k = nrow(schwein$map_loaded)-1, eig=T)
eigenvals(pcoa.bc.schwein)/sum(eigenvals(pcoa.bc.schwein)) # 17.8, 9.9
schwein$map_loaded$Axis01 <- scores(pcoa.bc.schwein)[,1]
schwein$map_loaded$Axis02 <- scores(pcoa.bc.schwein)[,2]
# c) troglodytes
pcoa.bc.trog <- cmdscale(trog.bacteria.distance, k = nrow(trog$map_loaded)-1, eig=T)
eigenvals(pcoa.bc.trog)/sum(eigenvals(pcoa.bc.trog)) # 13.6, 8.5
trog$map_loaded$Axis01 <- scores(pcoa.bc.trog)[,1]
trog$map_loaded$Axis02 <- scores(pcoa.bc.trog)[,2]
# d) verus
pcoa.bc.verus <- cmdscale(verus.bacteria.distance, k = nrow(verus$map_loaded)-1, eig=T)
eigenvals(pcoa.bc.verus)/sum(eigenvals(pcoa.bc.verus)) # 12.7, 6.2
verus$map_loaded$Axis01 <- scores(pcoa.bc.verus)[,1]
verus$map_loaded$Axis02 <- scores(pcoa.bc.verus)[,2]
# Combine
elli$map_loaded$dataset <- "a"
schwein$map_loaded$dataset <- "b"
trog$map_loaded$dataset <- "c"
verus$map_loaded$dataset <- "d"
pcoa_bc_multi <- rbind(elli$map_loaded, schwein$map_loaded, trog$map_loaded, verus$map_loaded)
pcoa_bc_multi$dataset <- as.factor(pcoa_bc_multi$dataset)

#### __Part 2: Jaccard ####
# a) ellioti
pcoa.jac.elli <- cmdscale(elli.parasites.distance, k = nrow(elli$map_loaded)-1, eig=T)
eigenvals(pcoa.jac.elli)/sum(eigenvals(pcoa.jac.elli)) # 47.7, 20.6
elli$map_loaded$Axis01p <- scores(pcoa.jac.elli)[,1]
elli$map_loaded$Axis02p <- scores(pcoa.jac.elli)[,2]
# b) schweinfurthii
pcoa.jac.schwein <- cmdscale(schwein.parasites.distance, k = nrow(schwein$map_loaded)-1, eig=T)
eigenvals(pcoa.jac.schwein)/sum(eigenvals(pcoa.jac.schwein)) # 27.5, 20.7
schwein$map_loaded$Axis01p <- scores(pcoa.jac.schwein)[,1]
schwein$map_loaded$Axis02p <- scores(pcoa.jac.schwein)[,2]
# c) troglodytes
pcoa.jac.trog <- cmdscale(trog.parasites.distance, k = nrow(trog$map_loaded)-1, eig=T)
eigenvals(pcoa.jac.trog)/sum(eigenvals(pcoa.jac.trog)) # 32.6, 19.8
trog$map_loaded$Axis01p <- scores(pcoa.jac.trog)[,1]
trog$map_loaded$Axis02p <- scores(pcoa.jac.trog)[,2]
# d) verus
pcoa.jac.verus <- cmdscale(verus.parasites.distance, k = nrow(verus$map_loaded)-1, eig=T)
eigenvals(pcoa.jac.verus)/sum(eigenvals(pcoa.jac.verus)) # 30.1, 21.6
verus$map_loaded$Axis01p <- scores(pcoa.jac.verus)[,1]
verus$map_loaded$Axis02p <- scores(pcoa.jac.verus)[,2]
# Combine
pcoa_jac_multi <- rbind(elli$map_loaded, schwein$map_loaded, trog$map_loaded, verus$map_loaded)
pcoa_jac_multi$dataset <- as.factor(pcoa_jac_multi$dataset)

#### __Part 3: Combine ####
df_bc <- pcoa_bc_multi %>%
  select(subspecies, dataset, Axis01, Axis02, Site) %>%
  mutate(taxon = "e")
df_bc$sampleID <- rownames(pcoa_bc_multi)
df_jac <- pcoa_jac_multi %>% 
  select(subspecies, dataset, Axis01p, Axis02p, Site) %>%
  rename(Axis01 = Axis01p,
         Axis02 = Axis02p) %>%
  mutate(taxon = "f")
df_jac$sampleID <- rownames(pcoa_jac_multi)
df_comb <- rbind(df_bc, df_jac)
df_comb$taxon <- as.factor(df_comb$taxon)
levels(df_comb$dataset)
df_comb$dataset <- factor(df_comb$dataset,
                          levels = c("d","a","c","b"))
sample_newdiet <- select(diet_update, sampleID, Diet)
df_comb <- left_join(df_comb, sample_newdiet, by = "sampleID")
facet_names.7 <- c('a' = "N-C",
                   'b' = "East",
                   'c' = "Central",
                   'd' = "West",
                   'e' = "Prokaryotes",
                   'f' = "Parasites")
label.df.7 <- data.frame(taxon = c("e","e","e","e","f","f","f","f"),
                         dataset = c("a","b","c","d","a","b","c","d"),
                         x = c(-0.5,-0.5,-0.5,-0.5,-0.65,-0.65,-0.65,-0.65),
                         y = c(-0.4,-0.55, -0.45, -0.25, -0.4,-0.55, -0.45, -0.25),
                         label = c("26.6%,10.3%","17.8%,9.9%","13.6%,8.5%","12.7%,6.2%",
                                   "47.7%,20.6%","27.5%,20.7%","32.6%,19.8%","30.1%,21.6%"))
# Add site centroids and names
# Make 8 centroid dataframes from betadisper (one for each panel), then combine
c1 <- betadisper(elli.bacteria.distance, elli$map_loaded$Site)
c1 <- data.frame(Site = rownames(c1$centroids),
                 Axis01 = c1$centroids[,1],
                 Axis02 = c1$centroids[,2],
                 dataset = "a",
                 taxon = "e")
c2 <- betadisper(schwein.bacteria.distance, schwein$map_loaded$Site)
c2 <- data.frame(Site = rownames(c2$centroids),
                 Axis01 = c2$centroids[,1],
                 Axis02 = c2$centroids[,2],
                 dataset = "b",
                 taxon = "e")
c3 <- betadisper(trog.bacteria.distance, trog$map_loaded$Site)
c3 <- data.frame(Site = rownames(c3$centroids),
                 Axis01 = c3$centroids[,1],
                 Axis02 = c3$centroids[,2],
                 dataset = "c",
                 taxon = "e")
c4 <- betadisper(verus.bacteria.distance, verus$map_loaded$Site)
c4 <- data.frame(Site = rownames(c4$centroids),
                 Axis01 = c4$centroids[,1],
                 Axis02 = c4$centroids[,2],
                 dataset = "d",
                 taxon = "e")
c5 <- betadisper(elli.parasites.distance, elli$map_loaded$Site)
c5 <- data.frame(Site = rownames(c5$centroids),
                 Axis01 = c5$centroids[,1],
                 Axis02 = c5$centroids[,2],
                 dataset = "a",
                 taxon = "f")
c6 <- betadisper(schwein.parasites.distance, schwein$map_loaded$Site)
c6 <- data.frame(Site = rownames(c6$centroids),
                 Axis01 = c6$centroids[,1],
                 Axis02 = c6$centroids[,2],
                 dataset = "b",
                 taxon = "f")
c7 <- betadisper(trog.parasites.distance, trog$map_loaded$Site)
c7 <- data.frame(Site = rownames(c7$centroids),
                 Axis01 = c7$centroids[,1],
                 Axis02 = c7$centroids[,2],
                 dataset = "c",
                 taxon = "f")
c8 <- betadisper(verus.parasites.distance, verus$map_loaded$Site)
c8 <- data.frame(Site = rownames(c8$centroids),
                 Axis01 = c8$centroids[,1],
                 Axis02 = c8$centroids[,2],
                 dataset = "d",
                 taxon = "f")
centroids <- rbind(c1,c2,c3,c4,c5,c6,c7,c8)
site_diet <- diet_update[ !duplicated(diet_update$Site), ] %>%
  select(Site, Diet)
centroids <- left_join(centroids, site_diet, by = "Site")
centroids.l <- subset(centroids, taxon == "e")
centroids.l <- subset(centroids.l, dataset == "a" | dataset == "b" | dataset == "c")
# Show site names. To make pretty just do in PPT.
# library(ggrepel)
# geom_text_repel(data = centroids.l,
#                 aes(label = Site),
#                 size = 3, box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) +
df_comb$Diet <- factor(df_comb$Diet,
                    levels = c("Termites", "Honey", "None", "Algae, Termites", "Honey, Termites",
                               "Nuts", "Nuts, Honey, Termites", "Algae", "Algae, Honey, Termites"))
centroids$Diet <- factor(centroids$Diet,
                    levels = c("Termites", "Honey", "None", "Algae, Termites", "Honey, Termites",
                               "Nuts", "Nuts, Honey, Termites", "Algae", "Algae, Honey, Termites"))
micro.hulls <- ddply(df_comb, c("taxon","dataset","Diet"), find_hull)
pdf("forPPT7.pdf", width = 6.5, height = 5.5)
ggplot(df_comb, aes(Axis01, Axis02, colour = Diet, fill = Diet)) +
  geom_polygon(data = micro.hulls, alpha = 0.1, show.legend = F) +
  geom_point(size = 1, alpha = 0.3) +
  geom_point(data = centroids, size = 1.5, pch = 17, show.legend = F) +
  geom_text(data = label.df.7, aes(x = x, y = y, label = label), size = 3, 
            inherit.aes = F) +
  labs(x = "PC1", y = "PC2", colour = "Specialty Items") +
  facet_grid(dataset ~ taxon, labeller = as_labeller(facet_names.7), scales = "free") +
  guides(colour = guide_legend(override.aes = list(alpha=1, size = 2)),
         fill = FALSE) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(c(0,-2,0,0))), 
        axis.title.x = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 12))
dev.off()



#### __Check Seasonality
verus_augjan <- filter_data()



########################## _8. Human - Chimp Geography Climate Sex ############################
# Add humans
humans <- select(humans, -X)
names(humans) == names(fig5comb)
humans_km <- select(humans_km, -X)
names(humans_km) == names(fig5comb)
fig5combGC <- subset(fig5comb, dataset == "Geography" | dataset == "Climate")
fig5combGCall <- subset(fig5combGC, subspecies == "all")
humans_clim <- subset(humans, dataset == "Climate")
humans_clim_km <- rbind(humans_clim, humans_km)
fig8 <- rbind(fig5combGCall, humans_clim_km)
(max(subset(fig8, dataset == "Geography")$dist) - min(subset(fig8, dataset == "Geography")$dist))/2
(max(subset(fig8, dataset == "Climate")$dist) - min(subset(fig8, dataset == "Climate")$dist))/2

label.df.8 <- data.frame(subspecies = c("all","all","human","human"),
                         dataset = c("Geography","Climate","Geography","Climate"),
                         x = c(2640.904, 3.602599,
                               2640.904, 3.602599),
                         y = c(0.125,0.125,0.125,0.125),
                         label = c("n = 560, r = 0.68, p = 0.001",
                                   "n = 560, r = 0.40, p = 0.001",
                                   "n = 2296, r = 0.007, p = 0.05",
                                   "n = 2296, r = 0.007, p = 0.29"))
facet_names.8 <- c("all" = "Chimpanzees", "human" = "Humans", 
                   "Geography" = "Geography", "Climate" = "Climate")
levels(fig8$dataset)
ggplot(data = fig8, aes(dist, BC)) +
  geom_point(data = subset(fig8, subspecies != "human"),
             size = 0.5, alpha = 0.01) +
  geom_point(data = subset(fig8, subspecies == "human"),
             size = 0.1, alpha = 0.01) +
  geom_smooth(data = subset(fig8, dataset == "Geography" & subspecies == "all"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig8, dataset == "Climate" & subspecies == "all"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig8, subspecies == "human"),
              method = lm, formula = y ~ x, size = 0.5, linetype = "dashed") +
  geom_text(data = label.df.8, aes(x = x, y = y, label = label), size = 3) +
  labs(x = "Distance",
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x", labeller = as_labeller(facet_names.8)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))
dev.off()

# Update: Compare Effect of Sex in humans vs. chimps too
# Run section 6 first - Similar analysis as within vs between sites
# Within F, within M, F v M
# Now add sex and subspecies, matching to col1 and col2
sex_subspecies <- select(bac_euk$map_loaded, sampleID, Sex, subspecies)
bac_bray_df_long_sex <- inner_join(bac_bray_df_long, sex_subspecies, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long_sex <- inner_join(bac_bray_df_long_sex, sex_subspecies, 
                               by = c("variable" = "sampleID"))
# Make new column indicating if within or between regions
for (i in 1:nrow(bac_bray_df_long_sex)) {
  ifelse(bac_bray_df_long_sex$subspecies.x[i] == bac_bray_df_long_sex$subspecies.y[i],
         bac_bray_df_long_sex$region[i] <- "within",
         bac_bray_df_long_sex$region[i] <- "between")
}

# Run if you want to just use within region comparisons
# bac_bray_df_long_sex <- subset(bac_bray_df_long_sex, region == "within")

# Make new column indicating if comparison is FF, MM or FM
# Vectorize and use ifelse for efficiency
c <- character(nrow(bac_bray_df_long_sex))
for (i in 1:nrow(bac_bray_df_long_sex)) {
  ifelse(bac_bray_df_long_sex$Sex.x[i] != bac_bray_df_long_sex$Sex.y[i],
         c[i] <- "FM",
         ifelse(bac_bray_df_long_sex$Sex.x[i]=="F"&bac_bray_df_long_sex$Sex.y[i]=="F",
                c[i] <- "FF",
                c[i] <- "MM"))
}
bac_bray_df_long_sex$comparison <- c
bac_bray_df_long_sex$comparison <- as.factor(bac_bray_df_long_sex$comparison)
levels(bac_bray_df_long_sex$comparison)
# Sample sizes
table(bac_bray_df_long_sex$comparison)
# Stats (all highly significant)
m.sex <- lm(value ~ comparison, data = bac_bray_df_long_sex)
summary(m.sex)
m.sex <- aov(value ~ comparison, data = bac_bray_df_long_sex)
summary(m.sex)
TukeyHSD(m.sex)
hist(m.sex$residuals)
leveneTest(value ~ comparison, data = bac_bray_df_long_sex)
kruskal.test(value ~ comparison, data = bac_bray_df_long_sex)
library(PMCMR)
posthoc.kruskal.nemenyi.test(value ~ comparison, data = bac_bray_df_long_sex)
bac_bray_df_long_sex$region2 <- bac_bray_df_long_sex$subspecies.x

# Graph
sex_means <- ddply(bac_bray_df_long_sex, c("comparison"), summarise,
                      mean = mean(value, na.rm = T),
                      se = se(value))
sex_means
label.df.sex <- data.frame(comparison = c("FF", "FM", "MM"),
                          Value = c(1.05,1.05,1.05),
                          Sig = c("a","b","c"))
label.df.sex.n <- data.frame(comparison = c("FF", "FM", "MM"),
                             Value = c(0,0,0),
                             Sig = c("n = 32896","n = 77871","n = 45753"))
label.df.sex.m <- data.frame(comparison = c("FF", "FM", "MM"),
                             Value = sex_means$mean - 0.035,
                             Sig = c(round(sex_means$mean[1], digits = 3),
                                     round(sex_means$mean[2], digits = 3),
                                     round(sex_means$mean[3], digits = 3)))
pdf("sex_chimp.pdf", width = 6.5, 4)
ggplot(data = bac_bray_df_long_sex, aes(comparison, value)) +
  geom_boxplot(outlier.alpha = 0.1, outlier.size = 1) +
  geom_point(data = sex_means, aes(comparison, mean), size = 4, shape = 17, color = "red") +
  geom_text(data = label.df.sex, aes(x = comparison, y = Value, label = Sig), size = 5) +
  geom_text(data = label.df.sex.n, aes(x = comparison, y = Value, label = Sig)) +
  geom_text(data = label.df.sex.m, aes(x = comparison, y = Value, label = Sig)) +
  labs(x = "Sex Comparison",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0,1.05) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
dev.off()

# Combine
# Facet_grid can't handle continuous and categorical so change FF, FM, MM to continuous
# Then in Preview cover cover labels and relabel with FF, FM, MM
names(fig8)
names(bac_bray_df_long_sex)
bac_bray_df_long_sex$subspecies <- "all"
bac_bray_df_long_sex$dataset <- "Sex"
names(humans_sex)
humans_sex$subspecies <- "human"
humans_sex$dataset <- "Sex"
chimps_sex <- dplyr::select(bac_bray_df_long_sex, sampleID, variable, value, subspecies, dataset, comparison)
humans_sex <- dplyr::select(humans_sex, sampleID, variable, value, subspecies, dataset, comparison)
names(chimps_sex) <- names(fig8)
names(humans_sex) <- names(fig8)
levels(humans_sex$dist)
humans_sex$dist <- as.character(humans_sex$dist)
humans_sex[humans_sex == "FF"] <- 2.5
humans_sex[humans_sex == "FM"] <- 7.5
humans_sex[humans_sex == "MM"] <- 12.5
levels(chimps_sex$dist)
chimps_sex$dist <- as.character(chimps_sex$dist)
chimps_sex[chimps_sex == "FF"] <- 2.5
chimps_sex[chimps_sex == "FM"] <- 7.5
chimps_sex[chimps_sex == "MM"] <- 12.5
humans_sex$dist <- as.numeric(humans_sex$dist)
chimps_sex$dist <- as.numeric(chimps_sex$dist)
str(fig8)
str(chimps_sex)
str(humans_sex)
chimps_sex$subspecies <- as.factor(chimps_sex$subspecies)
chimps_sex$dataset <- as.factor(chimps_sex$dataset)
chimps_sex$variable <- as.factor(chimps_sex$variable)
humans_sex$subspecies <- as.factor(humans_sex$subspecies)
humans_sex$dataset <- as.factor(humans_sex$dataset)
humans_sex$variable <- as.factor(humans_sex$variable)
humans_sex$sampleID <- as.character(humans_sex$variable)
fig8wsex <- rbind(fig8, chimps_sex, humans_sex)
str(fig8wsex)
levels(fig8wsex$dataset)
label.df.8.sex1 <- data.frame(subspecies = c("all","all","all","human","human","human"),
                         dataset = c("Geography","Climate","Sex","Geography","Climate","Sex"),
                         x = c(2640.904, 3.602599,7.5,
                               2640.904, 3.602599,7.5),
                         y = c(0.125,0.125,0.125,0.125,0.125,0.125),
                         label = c("n = 560, r = 0.68, p = 0.001",
                                   "n = 560, r = 0.40, p = 0.001",
                                   "n = 560, R^2 = 0.0008, p < 0.001",
                                   "n = 2296, r = 0.007, p = 0.05",
                                   "n = 2296, r = 0.007, p = 0.29",
                                   "n = 2218, R^2 = 0.0003, p < 0.001"))
facet_names.8.sex <- c("all" = "Chimpanzees", "human" = "Humans", 
                   "Geography" = "a) Geography", "Climate" = "b) Climate", "Sex" = "c) Sex")
lab.sex <- data.frame(subspecies = c("human","human","human"),
                      dataset = c("Sex","Sex","Sex"),
                      x = c(5,10,15),
                      y = c(0,0,0),
                      label = c("  ","  ","  "))
sex_both <- subset(fig8wsex, dataset == "Sex")
sex_means.both <- ddply(sex_both, c("subspecies","dist"), summarise,
                        mean = mean(BC, na.rm = T),
                        se = se(BC),
                        sd = sd(BC))
sex_means.both$dataset <- "Sex"
sex_means.both$dataset <- as.factor(sex_means.both$dataset)
tiff(file = "PNGs/Figure5.tiff", width = 5.6, height = 3.4, unit = "in", res = 300)
ggplot(data = fig8wsex, aes(dist, BC)) +
  geom_point(data = subset(fig8wsex, subspecies != "human" & dataset != "Sex"),
             size = 0.5, alpha = 0.01) +
  geom_point(data = subset(fig8wsex, subspecies == "human" & dataset != "Sex"),
             size = 0.1, alpha = 0.01) +
  geom_smooth(data = subset(fig8wsex, dataset == "Geography" & subspecies == "all"),
              method = lm, formula = y ~ poly(x, 2), size = 0.5) +
  geom_smooth(data = subset(fig8wsex, dataset == "Climate" & subspecies == "all"),
              method = lm, formula = y ~ x, size = 0.5) +
  geom_smooth(data = subset(fig8wsex, subspecies == "human" & dataset == "Geography"),
              method = lm, formula = y ~ x, size = 0.5, linetype = "dashed") +
  geom_smooth(data = subset(fig8wsex, subspecies == "human" & dataset == "Climate"),
              method = lm, formula = y ~ x, size = 0.5, linetype = "dashed") +
  geom_jitter(data = subset(fig8wsex, dataset == "Sex"), size = 0.1, alpha = 0.01) +
  geom_point(data = sex_means.both, aes(dist, mean), size = 3, shape = 16, color = "blue") +
  geom_errorbar(data = sex_means.both,
                aes(x = dist, ymin = mean-sd, ymax = mean+sd),
                width = 2, size = 0.5, inherit.aes = FALSE, color = "blue") +
  geom_text(data = label.df.8.sex1, aes(x = x, y = y, label = label), size = 2.25) +
  geom_text(data = lab.sex, aes(x = x, y = y, label = label), size = 5) +
  labs(x = NULL,
       y = "Bray-Curtis Dissimilarity") +
  facet_grid(subspecies ~ dataset, scales = "free_x", labeller=as_labeller(facet_names.8.sex)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text.y  = element_text(size = 12),
        axis.text.x = element_text(size = 12, margin = margin(c(0,0,0,0))),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))
dev.off()



#################################### _9. Genetics #############################################
# Genetics Update
# Reviewer wanted us to do more with genetics
# Jack Lester made new genetic distance matrices at the individual level
library(readxl)
sp.distance <- read_xlsx("../Genetics/smouse and peakall distances.xlsx") %>%
  column_to_rownames(var = "Obs.")
sp.distance.mat <- as.matrix(sp.distance)
sum(rownames(sp.distance.mat) != rownames(bacteria.distance.mat))
sum(colnames(sp.distance.mat) != colnames(bacteria.distance.mat))
sp.distance.mat <- sp.distance.mat[rownames(bacteria.distance.mat),
                                   rownames(bacteria.distance.mat)]
sum(rownames(sp.distance.mat) != rownames(bacteria.distance.mat))
sum(colnames(sp.distance.mat) != colnames(bacteria.distance.mat))
sp.distance.mat[upper.tri(sp.distance.mat, diag = TRUE)] <- NA
sp.distance <- as.dist(sp.distance.mat)

set.seed(500)
mantel(sp.distance.mat, bacteria.distance.mat) # r = 0.32, p = 0.001
pdf(file = "Space_Bray.pdf", width = 6.68, height = 4.46)
qplot(sp.distance, bacteria.distance, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Genetic Distance (Smouse and Peakall)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))
dev.off()

# Quickly look at the other methods Jack tested to see if they are similar
kin.distance.mat <- read_xlsx("../Genetics/kinship coeficient distances.xlsx") %>%
  column_to_rownames(var = "Obs.") %>%
  as.matrix()
kin.distance.mat <- kin.distance.mat[rownames(bacteria.distance.mat),
                                   rownames(bacteria.distance.mat)]
sum(rownames(kin.distance.mat) != rownames(bacteria.distance.mat))
sum(colnames(kin.distance.mat) != colnames(bacteria.distance.mat))
kin.distance.mat[upper.tri(kin.distance.mat, diag = TRUE)] <- NA
kin.distance <- as.dist(kin.distance.mat)
set.seed(500)
mantel(kin.distance.mat, bacteria.distance.mat) # r = -0.31, p = 1
qplot(kin.distance, bacteria.distance, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Genetic Distance (Kinship Coefficient)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

eucl.distance.mat <- read_xlsx("../Genetics/euclidean distances.xlsx") %>%
  column_to_rownames(var = "Obs.") %>%
  as.matrix()
eucl.distance.mat <- eucl.distance.mat[rownames(bacteria.distance.mat),
                                     rownames(bacteria.distance.mat)]
sum(rownames(eucl.distance.mat) != rownames(bacteria.distance.mat))
sum(colnames(eucl.distance.mat) != colnames(bacteria.distance.mat))
eucl.distance.mat[upper.tri(eucl.distance.mat, diag = TRUE)] <- NA
eucl.distance <- as.dist(eucl.distance.mat)
set.seed(500)
mantel(eucl.distance.mat, bacteria.distance.mat) # r = 0.34, p = 0.001
qplot(eucl.distance, bacteria.distance, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Genetic Distance (Kinship Coefficient)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

chord.distance.mat <- read_xlsx("../Genetics/chord distances.xlsx") %>%
  column_to_rownames(var = "Obs.") %>%
  as.matrix()
chord.distance.mat <- chord.distance.mat[rownames(bacteria.distance.mat),
                                       rownames(bacteria.distance.mat)]
sum(rownames(chord.distance.mat) != rownames(bacteria.distance.mat))
sum(colnames(chord.distance.mat) != colnames(bacteria.distance.mat))
chord.distance.mat[upper.tri(chord.distance.mat, diag = TRUE)] <- NA
chord.distance <- as.dist(chord.distance.mat)
set.seed(500)
mantel(chord.distance.mat, bacteria.distance.mat) # r = 0.40, p = 0.001
qplot(chord.distance, bacteria.distance, geom = c("point","smooth"), alpha = I(0.1)) +
  labs(x = "Genetic Distance (Chord)",
       y = "Bray-Curtis Dissimilarity") +
  ylim(0, 1) +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.1,0.1,0.1,0.15),"cm"))

#### __ New Figure S3 ####
# Need to make a multipanel figure of prok and euk correlations with genetic distance among and within regions (similar to the climate figure)
# Make regional matrices
elli.sp.distance.mat <- 
  sp.distance.mat[rownames(sp.distance.mat) %in% rownames(elli.bacteria.distance.mat),
                  rownames(sp.distance.mat) %in% rownames(elli.bacteria.distance.mat)]
schwein.sp.distance.mat <- 
  sp.distance.mat[rownames(sp.distance.mat) %in% rownames(schwein.bacteria.distance.mat),
                  rownames(sp.distance.mat) %in% rownames(schwein.bacteria.distance.mat)]
trog.sp.distance.mat <- 
  sp.distance.mat[rownames(sp.distance.mat) %in% rownames(trog.bacteria.distance.mat),
                  rownames(sp.distance.mat) %in% rownames(trog.bacteria.distance.mat)]
verus.sp.distance.mat <- 
  sp.distance.mat[rownames(sp.distance.mat) %in% rownames(verus.bacteria.distance.mat),
                  rownames(sp.distance.mat) %in% rownames(verus.bacteria.distance.mat)]


set.seed(500)
mantel(sp.distance.mat, bacteria.distance.mat) # r = 0.32, p = 0.001
mantel(elli.sp.distance.mat, elli.bacteria.distance.mat) # r = 0.42, p = 0.001
mantel(schwein.sp.distance.mat, schwein.bacteria.distance.mat) # r = 0.28, p = 0.001
mantel(trog.sp.distance.mat, trog.bacteria.distance.mat) # r = 0.22, p = 0.001
mantel(verus.sp.distance.mat, verus.bacteria.distance.mat) # r = 0.00, p = 0.51

mantel(sp.distance.mat, parasites.distance) # r = 0.21, p = 0.001
mantel(elli.sp.distance.mat, elli.parasites.distance) # r = 0.39, p = 0.002
mantel(schwein.sp.distance.mat, schwein.parasites.distance) # r = 0.03, p = 0.24
mantel(trog.sp.distance.mat, trog.parasites.distance) # r = 0.15, p = 0.001
mantel(verus.sp.distance.mat, verus.parasites.distance) # r = 0.01, p = 0.37

mantel.partial(sp.distance.mat, bacteria.distance.mat, 
               geography.distance.mat) # r 0.11,p 0.001
mantel.partial(elli.sp.distance.mat, elli.bacteria.distance.mat, 
               elli.geography.distance.mat) # r = 0.20, p = 0.03
mantel.partial(schwein.sp.distance.mat, schwein.bacteria.distance.mat, 
               schwein.geography.distance.mat) # r = 0.18, p = 0.001
mantel.partial(trog.sp.distance.mat, trog.bacteria.distance.mat, 
               trog.geography.distance.mat) # r = 0.10, p = 0.05
mantel.partial(verus.sp.distance.mat, verus.bacteria.distance.mat, 
               verus.geography.distance.mat) # r = -0.04, p = 0.91

mantel.partial(sp.distance.mat, parasites.distance.mat, 
               geography.distance.mat) # r = 0.06, p = 0.003
mantel.partial(elli.sp.distance.mat, elli.parasites.distance.mat, 
               elli.geography.distance.mat) # r = 0.18, p = 0.05
mantel.partial(schwein.sp.distance.mat, schwein.parasites.distance.mat, 
               schwein.geography.distance.mat) # r = 0.005, p = 0.45
mantel.partial(trog.sp.distance.mat, trog.parasites.distance.mat, 
               trog.geography.distance.mat) # r = 0.07, p = 0.08
mantel.partial(verus.sp.distance.mat, verus.parasites.distance.mat, 
               verus.geography.distance.mat) # r = -0.005, p = 0.59

# Make dataframes
sp.distance.df <- as.data.frame(sp.distance.mat)
sp.distance.df$sampleID <- rownames(sp.distance.df)
sp.distance.df.long <- melt(sp.distance.df, id.vars = "sampleID")
sp.distance.df.long <- na.omit(sp.distance.df.long)
bacteria.distance.df.long$SP <- sp.distance.df.long$value
parasites.distance.df.long$SP <- sp.distance.df.long$value

elli.sp.distance.df <- as.data.frame(elli.sp.distance.mat)
elli.sp.distance.df$sampleID <- rownames(elli.sp.distance.df)
elli.sp.distance.df.long <- melt(elli.sp.distance.df, id.vars = "sampleID")
elli.sp.distance.df.long <- na.omit(elli.sp.distance.df.long)
elli.bacteria.distance.df.long$SP <- elli.sp.distance.df.long$value
elli.parasites.distance.df.long$SP <- elli.sp.distance.df.long$value

schwein.sp.distance.df <- as.data.frame(schwein.sp.distance.mat)
schwein.sp.distance.df$sampleID <- rownames(schwein.sp.distance.df)
schwein.sp.distance.df.long <- melt(schwein.sp.distance.df, id.vars = "sampleID")
schwein.sp.distance.df.long <- na.omit(schwein.sp.distance.df.long)
schwein.bacteria.distance.df.long$SP <- schwein.sp.distance.df.long$value
schwein.parasites.distance.df.long$SP <- schwein.sp.distance.df.long$value

trog.sp.distance.df <- as.data.frame(trog.sp.distance.mat)
trog.sp.distance.df$sampleID <- rownames(trog.sp.distance.df)
trog.sp.distance.df.long <- melt(trog.sp.distance.df, id.vars = "sampleID")
trog.sp.distance.df.long <- na.omit(trog.sp.distance.df.long)
trog.bacteria.distance.df.long$SP <- trog.sp.distance.df.long$value
trog.parasites.distance.df.long$SP <- trog.sp.distance.df.long$value

verus.sp.distance.df <- as.data.frame(verus.sp.distance.mat)
verus.sp.distance.df$sampleID <- rownames(verus.sp.distance.df)
verus.sp.distance.df.long <- melt(verus.sp.distance.df, id.vars = "sampleID")
verus.sp.distance.df.long <- na.omit(verus.sp.distance.df.long)
verus.bacteria.distance.df.long$SP <- verus.sp.distance.df.long$value
verus.parasites.distance.df.long$SP <- verus.sp.distance.df.long$value

spgen.bact <- rbind(bacteria.distance.df.long,
                    elli.bacteria.distance.df.long,
                    trog.bacteria.distance.df.long,
                    schwein.bacteria.distance.df.long,
                    verus.bacteria.distance.df.long) %>%
  mutate(Dataset = "a) Prokaryotes")

spgen.euk <- rbind(parasites.distance.df.long,
                    elli.parasites.distance.df.long,
                    trog.parasites.distance.df.long,
                    schwein.parasites.distance.df.long,
                    verus.parasites.distance.df.long) %>%
  mutate(Dataset = "b) Parasites")

names(spgen.bact)[3] <- "dis"
names(spgen.euk)[3] <- "dis"

spgen <- rbind(spgen.bact, spgen.euk) %>%
  mutate_if(is.character, as.factor)
(max(sp.distance.df.long$value) - min(sp.distance.df.long$value))/2
spgen_text <- data.frame(subspecies = c("all","all", "ellioti", "ellioti", "schweinfurthii","schweinfurthii","troglodytes","troglodytes","verus","verus"),
                         Dataset = c("a) Prokaryotes", "b) Parasites",
                                     "a) Prokaryotes", "b) Parasites",
                                     "a) Prokaryotes", "b) Parasites",
                                     "a) Prokaryotes", "b) Parasites",
                                     "a) Prokaryotes", "b) Parasites"),
                         x = c(4, 4, 4, 4, 4,
                               4, 4, 4, 4, 4),
                         y = c(0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08),
                         label = c("n = 560, r = 0.32, p = 0.001",
                                   "n = 560, r = 0.21, p = 0.001",
                                   "n = 28, r = 0.42, p = 0.001",
                                   "n = 28, r = 0.39, p = 0.002",
                                   "n = 134, r = 0.28, p = 0.001",
                                   "n = 134, r = 0.03, p = 0.24",
                                   "n = 86, r = 0.22, p = 0.001",
                                   "n = 86, r = 0.15, p = 0.001",
                                   "n = 312, r = 0.00, p = 0.51",
                                   "n = 312, r = 0.01, p = 0.37"))
spgen_facetnames <- c("a) Prokaryotes" = "a) Prokaryotes", 
                       "b) Parasites" = "b) Parasites",
                       "all" = "all", "schweinfurthii" = "East", 
                       "troglodytes" = "Central", "verus" = "West",
                       "ellioti" = "N-C")
levels(spgen$subspecies)
spgen$subspecies <- factor(spgen$subspecies,
                           levels = c("all","verus","ellioti","troglodytes","schweinfurthii"))

png(file = "PNGs/SupplementaryFigureS3_Genetics.png", 
    width = 6, height = 6, units = "in", res = 300)
ggplot(data = spgen, aes(SP, dis)) +
  geom_point(size = 0.5, alpha = 0.02) +
  geom_smooth(data = subset(spgen, subspecies == "all" | subspecies == "ellioti" |
                              subspecies == "troglodytes"),
              method = lm, se = T, size = 0.5) +
  geom_smooth(data = subset(spgen, 
                            subspecies == "schweinfurthii" & Dataset == "a) Prokaryotes"),
              method = lm, se = T, size = 0.5) +
  geom_text(data = spgen_text, aes(x = x, y = y, label = label), size = 2.5) +
  labs(x = "Genetic Distance",
       y = "Dissimilarity") +
  facet_grid(subspecies ~ Dataset, labeller = as_labeller(spgen_facetnames)) +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        strip.text = element_text(size = 10))
dev.off()



#################################### Statistics ###############################################
#### __Permanova ####
# More Library issues - unload then reload vegan and RVAideMemoire
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(vegan)
library(RVAideMemoire)

# Other Variables
adonis(genetic.distance ~ pcoa_gen_df$subspecies, permutations = 999)
pairwise.perm.manova(genetic.distance, fact = pcoa_gen_df$subspecies, nperm = 999)
anova(betadisper(genetic.distance, pcoa_gen_df$subspecies))

adonis(sl.diet.distance ~ sitemap$subspecies, permutations = 999) # p = 0.07, F = 1.6, R2 = 0.16
pairwise.perm.manova(sl.diet.distance, fact = sitemap$subspecies, nperm = 999) # NS
anova(betadisper(sl.diet.distance, sitemap$subspecies)) # NS

adonis(sl.climate.distance ~ sitemap$subspecies, permutations = 999) # p=0.001,F = 4.5, R2 = 0.35
pairwise.perm.manova(sl.climate.distance, fact = sitemap$subspecies, nperm = 999)
anova(betadisper(sl.climate.distance, sitemap$subspecies)) # NS

set.seed(500)
adonis(plantrel ~ plantmeta$subspecies, permutations = 999)
pairwise.perm.manova(plantbray.dist, fact = plantmeta$subspecies, nperm = 999)
anova(betadisper(plantbray.dist, plantmeta$subspecies))
pcoa.veg <- cmdscale(plantbray.dist, k = nrow(plantmeta)-1, eig=T)

# Microbes
set.seed(500)
adonis(bacteria.distance ~ bac_euk$map_loaded$habitat + bac_euk$map_loaded$Diet + bac_euk$map_loaded$Site + bac_euk$map_loaded$Sex + bac_euk$map_loaded$Site:bac_euk$map_loaded$Sex, strata = bac_euk$map_loaded$subspecies)
set.seed(500)
adonis(bacteria.distance ~ bac_euk$map_loaded$habitat + bac_euk$map_loaded$DietD + bac_euk$map_loaded$Site + bac_euk$map_loaded$Sex + bac_euk$map_loaded$Site:bac_euk$map_loaded$Sex, strata = bac_euk$map_loaded$subspecies)

set.seed(500)
adonis(parasites.distance ~ bac_euk$map_loaded$habitat + bac_euk$map_loaded$Diet + bac_euk$map_loaded$Site + bac_euk$map_loaded$Sex + bac_euk$map_loaded$Site:bac_euk$map_loaded$Sex, strata = bac_euk$map_loaded$subspecies)
set.seed(500)
adonis(parasites.distance ~ bac_euk$map_loaded$habitat + bac_euk$map_loaded$DietD + bac_euk$map_loaded$Site + bac_euk$map_loaded$Sex + bac_euk$map_loaded$Site:bac_euk$map_loaded$Sex, strata = bac_euk$map_loaded$subspecies)



#### __GDM ####
library(dplyr)
library(TSdist)
library(raster) # Note this will mess up the select function in dplyr
library(gdm)
# Make matrices for the veg subset data
geography.distance.v.mat <- as.matrix(geog.distance.v)
climate.v <- dplyr::select(bac_euk_veg$map_loaded, 
                           AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
climate.v <- as.data.frame(scale(climate.v))
climate.v.distance <- dist(climate.v, method = "euclidean")
elli.climate.v <- dplyr::select(elli_veg$map_loaded, 
                                AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
elli.climate.v <- as.data.frame(scale(elli.climate.v))
elli.climate.v.distance <- dist(elli.climate.v, method = "euclidean")
schwein.climate.v <- dplyr::select(schwein_veg$map_loaded,
                                   AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
schwein.climate.v <- as.data.frame(scale(schwein.climate.v))
schwein.climate.v.distance <- dist(schwein.climate.v, method = "euclidean")
trog.climate.v <- dplyr::select(trog_veg$map_loaded, 
                                AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
trog.climate.v <- as.data.frame(scale(trog.climate.v))
trog.climate.v.distance <- dist(trog.climate.v, method = "euclidean")
verus.climate.v <- dplyr::select(verus_veg$map_loaded, 
                                 AnTemp, AnPercip, Seasnlty, PrcpSeasnlty)
verus.climate.v <- as.data.frame(scale(verus.climate.v))
verus.climate.v.distance <- dist(verus.climate.v, method = "euclidean")
climate.v.distance.mat <- as.matrix(climate.v.distance)
elli.climate.v.distance.mat <- as.matrix(elli.climate.v.distance)
schwein.climate.v.distance.mat <- as.matrix(schwein.climate.v.distance)
trog.climate.v.distance.mat <- as.matrix(trog.climate.v.distance)
verus.climate.v.distance.mat <- as.matrix(verus.climate.v.distance)
# Diet Jaccard. If warning given, must change NA to 0 or 1 (for sites with all zeroes)
diet.v <- dplyr::select(bac_euk_veg$map_loaded, 
                        algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                        termites, tubers, water)
diet.v.distance <- vegdist(diet.v, method = "jaccard")
diet.v.distance <- dist.zeroes(diet.v, diet.v.distance)
elli.diet.v <- dplyr::select(elli_veg$map_loaded, 
                             algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                             termites, tubers, water)
elli.diet.v.distance <- vegdist(elli.diet.v, method = "jaccard")
schwein.diet.v <- dplyr::select(schwein_veg$map_loaded,
                                algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                                termites, tubers, water)
schwein.diet.v.distance <- vegdist(schwein.diet.v, method = "jaccard")
schwein.diet.v.distance <- dist.zeroes(schwein.diet.v, schwein.diet.v.distance)
trog.diet.v <- dplyr::select(trog_veg$map_loaded, 
                             algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                             termites, tubers, water)
trog.diet.v.distance <- vegdist(trog.diet.v, method = "jaccard")
verus.diet.v <- dplyr::select(verus_veg$map_loaded, 
                              algae, ants, fruit, honey, marrow, meat, nuts, 'palm heart', 
                              termites, tubers, water)
verus.diet.v.distance <- vegdist(verus.diet.v, method = "jaccard")
verus.diet.v.distance <- dist.zeroes(verus.diet.v, verus.diet.v.distance)
diet.v.distance.mat <- as.matrix(diet.v.distance)
elli.diet.v.distance.mat <- as.matrix(elli.diet.v.distance)
schwein.diet.v.distance.mat <- as.matrix(schwein.diet.v.distance)
trog.diet.v.distance.mat <- as.matrix(trog.diet.v.distance)
verus.diet.v.distance.mat <- as.matrix(verus.diet.v.distance)



#### ____Prokaryotes ####
gdm.sampID <- rownames(bacteria.distance.v.mat)
bacteria.distance.v.mat.gdm <- cbind(gdm.sampID, bacteria.distance.v.mat)
geography.distance.v.mat.gdm <- cbind(gdm.sampID, geography.distance.v.mat)
climate.distance.v.mat.gdm <- cbind(gdm.sampID, climate.v.distance.mat)
veg.distance.mat.gdm <- cbind(gdm.sampID, veg.distance.mat)
diet.distance.v.mat.gdm <- cbind(gdm.sampID, diet.v.distance.mat)
bac_euk_veg$map_loaded$gdm.sampID <- bac_euk_veg$map_loaded$sampleID
gdm.data <- dplyr::select(bac_euk_veg$map_loaded, gdm.sampID, lat, long)
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = gdm.data,
                           distPreds = list(climate.distance.v.mat.gdm,
                                            veg.distance.mat.gdm,
                                            diet.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE)
summary(gdm.1)
str(gdm.1)
length(gdm.1$predictors)
plot(gdm.1, plot.layout = c(2,3))
gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
par(mfrow = c(1,4))
plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_1"], gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_2"], gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_3"], gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
max(gdm.1.splineDat$y[,"matrix_3"])
gdm.1.pred <- predict(gdm.1, gdm.bray)
head(gdm.1.pred)
par(mfrow = c(1,1))
plot(gdm.bray$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# schweinfurthii
schwein.gdm.sampID <- rownames(schwein.bacteria.distance.v.mat)
schwein.bacteria.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.bacteria.distance.v.mat)
schwein.climate.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.climate.v.distance.mat)
schwein.veg.distance.mat.gdm <- cbind(schwein.gdm.sampID, schwein.veg.distance.mat)
schwein.diet.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.diet.v.distance.mat)
schwein_veg$map_loaded$schwein.gdm.sampID <- schwein_veg$map_loaded$sampleID
schwein.gdm.data <- dplyr::select(schwein_veg$map_loaded, schwein.gdm.sampID, lat, long)
schwein.gdm.bray <- formatsitepair(bioData = schwein.bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "schwein.gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = schwein.gdm.data,
                           distPreds = list(schwein.climate.distance.v.mat.gdm,
                                            schwein.veg.distance.mat.gdm,
                                            schwein.diet.distance.v.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.bray, geo = TRUE)
summary(schwein.gdm.1)
str(schwein.gdm.1)
length(schwein.gdm.1$predictors)
plot(schwein.gdm.1, plot.layout = c(2,3))
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
str(schwein.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(schwein.gdm.1.splineDat$x[,"Geographic"], schwein.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_1"], schwein.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_2"], schwein.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_3"], schwein.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
max(schwein.gdm.1.splineDat$y[,"matrix_2"])
max(schwein.gdm.1.splineDat$y[,"matrix_3"])
schwein.gdm.1.pred <- predict(schwein.gdm.1, schwein.gdm.bray)
head(schwein.gdm.1.pred)
par(mfrow = c(1,1))
plot(schwein.gdm.bray$distance, schwein.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# troglodytes
trog.gdm.sampID <- rownames(trog.bacteria.distance.v.mat)
trog.bacteria.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.bacteria.distance.v.mat)
trog.climate.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.climate.v.distance.mat)
trog.veg.distance.mat.gdm <- cbind(trog.gdm.sampID, trog.veg.distance.mat)
trog.diet.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.diet.v.distance.mat)
trog_veg$map_loaded$trog.gdm.sampID <- trog_veg$map_loaded$sampleID
trog.gdm.data <- dplyr::select(trog_veg$map_loaded, trog.gdm.sampID, lat, long)
trog.gdm.bray <- formatsitepair(bioData = trog.bacteria.distance.v.mat.gdm,
                                   bioFormat = 3,
                                   siteColumn = "trog.gdm.sampID",
                                   XColumn = "long",
                                   YColumn = "lat",
                                   predData = trog.gdm.data,
                                   distPreds = list(trog.climate.distance.v.mat.gdm,
                                                    trog.veg.distance.mat.gdm,
                                                    trog.diet.distance.v.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.bray, geo = TRUE)
summary(trog.gdm.1)
str(trog.gdm.1)
length(trog.gdm.1$predictors)
plot(trog.gdm.1, plot.layout = c(2,3))
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
str(trog.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(trog.gdm.1.splineDat$x[,"Geographic"], trog.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_1"], trog.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_2"], trog.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_3"], trog.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
max(trog.gdm.1.splineDat$y[,"matrix_2"])
max(trog.gdm.1.splineDat$y[,"matrix_3"])
trog.gdm.1.pred <- predict(trog.gdm.1, trog.gdm.bray)
head(trog.gdm.1.pred)
par(mfrow = c(1,1))
plot(trog.gdm.bray$distance, trog.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# verus
verus.gdm.sampID <- rownames(verus.bacteria.distance.v.mat)
verus.bacteria.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.bacteria.distance.v.mat)
verus.climate.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.climate.v.distance.mat)
verus.veg.distance.mat.gdm <- cbind(verus.gdm.sampID, verus.veg.distance.mat)
verus.diet.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.diet.v.distance.mat)
verus_veg$map_loaded$verus.gdm.sampID <- verus_veg$map_loaded$sampleID
verus.gdm.data <- dplyr::select(verus_veg$map_loaded, verus.gdm.sampID, lat, long)
verus.gdm.bray <- formatsitepair(bioData = verus.bacteria.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "verus.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = verus.gdm.data,
                                distPreds = list(verus.climate.distance.v.mat.gdm,
                                                 verus.veg.distance.mat.gdm,
                                                 verus.diet.distance.v.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.bray, geo = TRUE)
summary(verus.gdm.1)
str(verus.gdm.1)
length(verus.gdm.1$predictors)
plot(verus.gdm.1, plot.layout = c(2,3))
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
str(verus.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(verus.gdm.1.splineDat$x[,"Geographic"], verus.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_1"], verus.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_2"], verus.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_3"], verus.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
max(verus.gdm.1.splineDat$y[,"matrix_2"])
max(verus.gdm.1.splineDat$y[,"matrix_3"])
verus.gdm.1.pred <- predict(verus.gdm.1, verus.gdm.bray)
head(verus.gdm.1.pred)
par(mfrow = c(1,1))
plot(verus.gdm.bray$distance, verus.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



#### ____Parasites ####
# Parasites matrices made in Habitat section _8.

# Whole dataset
gdm.sampID <- rownames(parasites.v.distance.mat)
parasites.distance.v.mat.gdm <- cbind(gdm.sampID, parasites.v.distance.mat)
geography.distance.v.mat.gdm <- cbind(gdm.sampID, geography.distance.v.mat)
climate.distance.v.mat.gdm <- cbind(gdm.sampID, climate.v.distance.mat)
veg.distance.mat.gdm <- cbind(gdm.sampID, veg.distance.mat)
diet.distance.v.mat.gdm <- cbind(gdm.sampID, diet.v.distance.mat)
bac_euk_veg$map_loaded$gdm.sampID <- bac_euk_veg$map_loaded$sampleID
gdm.data <- dplyr::select(bac_euk_veg$map_loaded, gdm.sampID, lat, long)
gdm.jac <- formatsitepair(bioData = parasites.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = gdm.data,
                           distPreds = list(climate.distance.v.mat.gdm,
                                            veg.distance.mat.gdm,
                                            diet.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE)
summary(gdm.1)
str(gdm.1)
length(gdm.1$predictors)
plot(gdm.1, plot.layout = c(2,3))
gdm.1.splineDat <- isplineExtract(gdm.1)
str(gdm.1.splineDat)
par(mfrow = c(1,4))
plot(gdm.1.splineDat$x[,"Geographic"], gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_1"], gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_2"], gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(gdm.1.splineDat$x[,"matrix_3"], gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
max(gdm.1.splineDat$y[,"matrix_3"])
gdm.1.pred <- predict(gdm.1, gdm.jac)
head(gdm.1.pred)
par(mfrow = c(1,1))
plot(gdm.jac$distance, gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# schweinfurthii
schwein.gdm.sampID <- rownames(schwein.parasites.v.distance.mat)
schwein.parasites.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.parasites.v.distance.mat)
schwein.climate.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.climate.v.distance.mat)
schwein.veg.distance.mat.gdm <- cbind(schwein.gdm.sampID, schwein.veg.distance.mat)
schwein.diet.distance.v.mat.gdm <- cbind(schwein.gdm.sampID, schwein.diet.v.distance.mat)
schwein_veg$map_loaded$schwein.gdm.sampID <- schwein_veg$map_loaded$sampleID
schwein.gdm.data <- dplyr::select(schwein_veg$map_loaded, schwein.gdm.sampID, lat, long)
schwein.gdm.jac <- formatsitepair(bioData = schwein.parasites.distance.v.mat.gdm,
                                   bioFormat = 3,
                                   siteColumn = "schwein.gdm.sampID",
                                   XColumn = "long",
                                   YColumn = "lat",
                                   predData = schwein.gdm.data,
                                   distPreds = list(schwein.climate.distance.v.mat.gdm,
                                                    schwein.veg.distance.mat.gdm,
                                                    schwein.diet.distance.v.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.jac, geo = TRUE)
summary(schwein.gdm.1)
str(schwein.gdm.1)
length(schwein.gdm.1$predictors)
plot(schwein.gdm.1, plot.layout = c(2,3))
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
str(schwein.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(schwein.gdm.1.splineDat$x[,"Geographic"], schwein.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_1"], schwein.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_2"], schwein.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(schwein.gdm.1.splineDat$x[,"matrix_3"], schwein.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
max(schwein.gdm.1.splineDat$y[,"matrix_2"])
max(schwein.gdm.1.splineDat$y[,"matrix_3"])
schwein.gdm.1.pred <- predict(schwein.gdm.1, schwein.gdm.jac)
head(schwein.gdm.1.pred)
par(mfrow = c(1,1))
plot(schwein.gdm.jac$distance, schwein.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# troglodytes
trog.gdm.sampID <- rownames(trog.parasites.v.distance.mat)
trog.parasites.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.parasites.v.distance.mat)
trog.climate.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.climate.v.distance.mat)
trog.veg.distance.mat.gdm <- cbind(trog.gdm.sampID, trog.veg.distance.mat)
trog.diet.distance.v.mat.gdm <- cbind(trog.gdm.sampID, trog.diet.v.distance.mat)
trog_veg$map_loaded$trog.gdm.sampID <- trog_veg$map_loaded$sampleID
trog.gdm.data <- dplyr::select(trog_veg$map_loaded, trog.gdm.sampID, lat, long)
trog.gdm.jac <- formatsitepair(bioData = trog.parasites.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "trog.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = trog.gdm.data,
                                distPreds = list(trog.climate.distance.v.mat.gdm,
                                                 trog.veg.distance.mat.gdm,
                                                 trog.diet.distance.v.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.jac, geo = TRUE)
summary(trog.gdm.1)
str(trog.gdm.1)
length(trog.gdm.1$predictors)
plot(trog.gdm.1, plot.layout = c(2,3))
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
str(trog.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(trog.gdm.1.splineDat$x[,"Geographic"], trog.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_1"], trog.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_2"], trog.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(trog.gdm.1.splineDat$x[,"matrix_3"], trog.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
max(trog.gdm.1.splineDat$y[,"matrix_2"])
max(trog.gdm.1.splineDat$y[,"matrix_3"])
trog.gdm.1.pred <- predict(trog.gdm.1, trog.gdm.jac)
head(trog.gdm.1.pred)
par(mfrow = c(1,1))
plot(trog.gdm.jac$distance, trog.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

# verus
verus.gdm.sampID <- rownames(verus.parasites.v.distance.mat)
verus.parasites.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.parasites.v.distance.mat)
verus.climate.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.climate.v.distance.mat)
verus.veg.distance.mat.gdm <- cbind(verus.gdm.sampID, verus.veg.distance.mat)
verus.diet.distance.v.mat.gdm <- cbind(verus.gdm.sampID, verus.diet.v.distance.mat)
verus_veg$map_loaded$verus.gdm.sampID <- verus_veg$map_loaded$sampleID
verus.gdm.data <- dplyr::select(verus_veg$map_loaded, verus.gdm.sampID, lat, long)
verus.gdm.jac <- formatsitepair(bioData = verus.parasites.distance.v.mat.gdm,
                                 bioFormat = 3,
                                 siteColumn = "verus.gdm.sampID",
                                 XColumn = "long",
                                 YColumn = "lat",
                                 predData = verus.gdm.data,
                                 distPreds = list(verus.climate.distance.v.mat.gdm,
                                                  verus.veg.distance.mat.gdm,
                                                  verus.diet.distance.v.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.jac, geo = TRUE)
summary(verus.gdm.1)
str(verus.gdm.1)
length(verus.gdm.1$predictors)
plot(verus.gdm.1, plot.layout = c(2,3))
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
str(verus.gdm.1.splineDat)
par(mfrow = c(1,4))
plot(verus.gdm.1.splineDat$x[,"Geographic"], verus.gdm.1.splineDat$y[,"Geographic"], lwd=3,
     type="l", xlab="Geographic distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_1"], verus.gdm.1.splineDat$y[,"matrix_1"], lwd=3,
     type="l", xlab="Climate distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_2"], verus.gdm.1.splineDat$y[,"matrix_2"], lwd=3,
     type="l", xlab="Veg distance", ylab="Partial ecological distance")
plot(verus.gdm.1.splineDat$x[,"matrix_3"], verus.gdm.1.splineDat$y[,"matrix_3"], lwd=3,
     type="l", xlab="Diet distance", ylab="Partial ecological distance")
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
max(verus.gdm.1.splineDat$y[,"matrix_2"])
max(verus.gdm.1.splineDat$y[,"matrix_3"])
verus.gdm.1.pred <- predict(verus.gdm.1, verus.gdm.jac)
head(verus.gdm.1.pred)
par(mfrow = c(1,1))
plot(verus.gdm.jac$distance, verus.gdm.1.pred, xlab="Observed dissimilarity",
     ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))



#### ____No Clim ####
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = gdm.data,
                           distPreds = list(veg.distance.mat.gdm,
                                            diet.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE)
gdm.1.splineDat <- isplineExtract(gdm.1)
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
round(gdm.1$explained, digits = 2)

schwein.gdm.bray <- formatsitepair(bioData = schwein.bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "schwein.gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = schwein.gdm.data,
                           distPreds = list(schwein.veg.distance.mat.gdm,
                                            schwein.diet.distance.v.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.bray, geo = TRUE)
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
max(schwein.gdm.1.splineDat$y[,"matrix_2"])
round(schwein.gdm.1$explained, digits = 2)

trog.gdm.bray <- formatsitepair(bioData = trog.bacteria.distance.v.mat.gdm,
                                   bioFormat = 3,
                                   siteColumn = "trog.gdm.sampID",
                                   XColumn = "long",
                                   YColumn = "lat",
                                   predData = trog.gdm.data,
                                   distPreds = list(trog.veg.distance.mat.gdm,
                                                    trog.diet.distance.v.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.bray, geo = TRUE)
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
max(trog.gdm.1.splineDat$y[,"matrix_2"])
round(trog.gdm.1$explained, digits = 2)

verus.gdm.bray <- formatsitepair(bioData = verus.bacteria.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "verus.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = verus.gdm.data,
                                distPreds = list(verus.veg.distance.mat.gdm,
                                                 verus.diet.distance.v.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.bray, geo = TRUE)
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
max(verus.gdm.1.splineDat$y[,"matrix_2"])
round(verus.gdm.1$explained, digits = 2)



gdm.jac <- formatsitepair(bioData = parasites.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = gdm.data,
                           distPreds = list(veg.distance.mat.gdm,
                                            diet.distance.v.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE)
gdm.1.splineDat <- isplineExtract(gdm.1)
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
max(gdm.1.splineDat$y[,"matrix_2"])
round(gdm.1$explained, digits = 2)

schwein.gdm.jac <- formatsitepair(bioData = schwein.parasites.distance.v.mat.gdm,
                                   bioFormat = 3,
                                   siteColumn = "schwein.gdm.sampID",
                                   XColumn = "long",
                                   YColumn = "lat",
                                   predData = schwein.gdm.data,
                                   distPreds = list(schwein.veg.distance.mat.gdm,
                                                    schwein.diet.distance.v.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.jac, geo = TRUE)
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
max(schwein.gdm.1.splineDat$y[,"matrix_2"])
round(schwein.gdm.1$explained, digits = 2)

trog.gdm.jac <- formatsitepair(bioData = trog.parasites.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "trog.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = trog.gdm.data,
                                distPreds = list(trog.veg.distance.mat.gdm,
                                                 trog.diet.distance.v.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.jac, geo = TRUE)
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
max(trog.gdm.1.splineDat$y[,"matrix_2"])
round(trog.gdm.1$explained, digits = 2)

verus.gdm.jac <- formatsitepair(bioData = verus.parasites.distance.v.mat.gdm,
                                 bioFormat = 3,
                                 siteColumn = "verus.gdm.sampID",
                                 XColumn = "long",
                                 YColumn = "lat",
                                 predData = verus.gdm.data,
                                 distPreds = list(verus.veg.distance.mat.gdm,
                                                  verus.diet.distance.v.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.jac, geo = TRUE)
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
max(verus.gdm.1.splineDat$y[,"matrix_2"])
round(verus.gdm.1$explained, digits = 2)
detach("package:gdm", unload=TRUE)
detach("package:raster", unload=TRUE)

#### ____No Diet ####
gdm.bray <- formatsitepair(bioData = bacteria.distance.v.mat.gdm,
                           bioFormat = 3,
                           siteColumn = "gdm.sampID",
                           XColumn = "long",
                           YColumn = "lat",
                           predData = gdm.data,
                           distPreds = list(veg.distance.mat.gdm))
gdm.1 <- gdm(gdm.bray, geo = TRUE)
gdm.1.splineDat <- isplineExtract(gdm.1)
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
round(gdm.1$explained, digits = 2)

schwein.gdm.bray <- formatsitepair(bioData = schwein.bacteria.distance.v.mat.gdm,
                                   bioFormat = 3,
                                   siteColumn = "schwein.gdm.sampID",
                                   XColumn = "long",
                                   YColumn = "lat",
                                   predData = schwein.gdm.data,
                                   distPreds = list(schwein.veg.distance.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.bray, geo = TRUE)
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
round(schwein.gdm.1$explained, digits = 2)

trog.gdm.bray <- formatsitepair(bioData = trog.bacteria.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "trog.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = trog.gdm.data,
                                distPreds = list(trog.veg.distance.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.bray, geo = TRUE)
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
round(trog.gdm.1$explained, digits = 2)

verus.gdm.bray <- formatsitepair(bioData = verus.bacteria.distance.v.mat.gdm,
                                 bioFormat = 3,
                                 siteColumn = "verus.gdm.sampID",
                                 XColumn = "long",
                                 YColumn = "lat",
                                 predData = verus.gdm.data,
                                 distPreds = list(verus.veg.distance.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.bray, geo = TRUE)
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
round(verus.gdm.1$explained, digits = 2)



gdm.jac <- formatsitepair(bioData = parasites.distance.v.mat.gdm,
                          bioFormat = 3,
                          siteColumn = "gdm.sampID",
                          XColumn = "long",
                          YColumn = "lat",
                          predData = gdm.data,
                          distPreds = list(veg.distance.mat.gdm))
gdm.1 <- gdm(gdm.jac, geo = TRUE)
gdm.1.splineDat <- isplineExtract(gdm.1)
max(gdm.1.splineDat$y[,"Geographic"])
max(gdm.1.splineDat$y[,"matrix_1"])
round(gdm.1$explained, digits = 2)

schwein.gdm.jac <- formatsitepair(bioData = schwein.parasites.distance.v.mat.gdm,
                                  bioFormat = 3,
                                  siteColumn = "schwein.gdm.sampID",
                                  XColumn = "long",
                                  YColumn = "lat",
                                  predData = schwein.gdm.data,
                                  distPreds = list(schwein.veg.distance.mat.gdm))
schwein.gdm.1 <- gdm(schwein.gdm.jac, geo = TRUE)
schwein.gdm.1.splineDat <- isplineExtract(schwein.gdm.1)
max(schwein.gdm.1.splineDat$y[,"Geographic"])
max(schwein.gdm.1.splineDat$y[,"matrix_1"])
round(schwein.gdm.1$explained, digits = 2)

trog.gdm.jac <- formatsitepair(bioData = trog.parasites.distance.v.mat.gdm,
                               bioFormat = 3,
                               siteColumn = "trog.gdm.sampID",
                               XColumn = "long",
                               YColumn = "lat",
                               predData = trog.gdm.data,
                               distPreds = list(trog.veg.distance.mat.gdm))
trog.gdm.1 <- gdm(trog.gdm.jac, geo = TRUE)
trog.gdm.1.splineDat <- isplineExtract(trog.gdm.1)
max(trog.gdm.1.splineDat$y[,"Geographic"])
max(trog.gdm.1.splineDat$y[,"matrix_1"])
round(trog.gdm.1$explained, digits = 2)

verus.gdm.jac <- formatsitepair(bioData = verus.parasites.distance.v.mat.gdm,
                                bioFormat = 3,
                                siteColumn = "verus.gdm.sampID",
                                XColumn = "long",
                                YColumn = "lat",
                                predData = verus.gdm.data,
                                distPreds = list(verus.veg.distance.mat.gdm))
verus.gdm.1 <- gdm(verus.gdm.jac, geo = TRUE)
verus.gdm.1.splineDat <- isplineExtract(verus.gdm.1)
max(verus.gdm.1.splineDat$y[,"Geographic"])
max(verus.gdm.1.splineDat$y[,"matrix_1"])
round(verus.gdm.1$explained, digits = 2)
detach("package:gdm", unload=TRUE)
detach("package:raster", unload=TRUE)



#################################### End Script ###############################################