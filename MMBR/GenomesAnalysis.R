# Methylotrophic Genome Analysis
# IMG search of methylotrophic genes
# Look at number of genomes, species, genera
# Look at how many have each gene in pathway

#### Setup ####
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(scales)
`%notin%` <- Negate(`%in%`)

setwd("~/Desktop/Review/")



#### Copy number check ####
# Import KO function profile table
# Filter out the dataset name (# of genomes with the KO)
# Filter out the KOs of interest
# Summarize by KO and calculate max and min copy number
# Update with only high quality genomes (5/11/22)
ko_table <- read.delim("mcrA_high_KO_profile/profile.txt",
                       header = F) %>%
  filter(V1 != "mcrA_high") %>%
  filter(V2 == "KO:K00399" | V2 == "KO:K00194" | V2 == "KO:K00440" | V2 == "KO:K14080" |
           V2 == "KO:K14084" | V2 == "KO:K16179" | V2 == "KO:K16177" | V2 == "KO:K16954") %>%
  group_by(V2) %>%
  summarise(min_copy_num = min(V3),
            max_copy_num = max(V3))



#### Bacterial/Archaeal Count ####
# Do for each gene in the metagenome survey graph
# All: mcrA, K00399
# Ac: cdhD, K00194
# H2: mtrA, K00577 (old)
# H2: frhA, K00440 (updated)
# M: mtaA, K14080
# TMA: mttC, K14084
# DMA: mtbC, K16179
# MMA: mtmC, K16177
# DMS/MeSH/MMPA: mtsA, K16954
mcrA <- read.delim("exportdata_mcrA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mcrA_dom <- mcrA %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "All")
mcrA_gen <- mcrA %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "All")
mcrA_sp <- mcrA %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "All")
mcrA_dom_gen <- mcrA %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "All")
  
cdhD <- read.delim("exportdata_cdhD_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
cdhD_dom <- cdhD %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Acetate")
cdhD_gen <- cdhD %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Acetate")
cdhD_sp <- cdhD %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Acetate")
cdhD_dom_gen <- cdhD %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Acetate")

mtrA <- read.delim("exportdata_mtrA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtrA_dom <- mtrA %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
mtrA_gen <- mtrA %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
mtrA_sp <- mtrA %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
mtrA_dom_gen <- mtrA %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")

frhA <- read.delim("exportdata_frhA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
frhA_dom <- frhA %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
frhA_gen <- frhA %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
frhA_sp <- frhA %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")
frhA_dom_gen <- frhA %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "H2/CO2")

mtaA <- read.delim("exportdata_mtaA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtaA_dom <- mtaA %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
mtaA_gen <- mtaA %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
mtaA_sp <- mtaA %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
mtaA_dom_gen <- mtaA %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")

mttC <- read.delim("exportdata_mttC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mttC_dom <- mttC %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "TMA")
mttC_gen <- mttC %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "TMA")
mttC_sp <- mttC %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "TMA")
mttC_dom_gen <- mttC %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "TMA")

mtbC <- read.delim("exportdata_mtbC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtbC_dom <- mtbC %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMA")
mtbC_gen <- mtbC %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMA")
mtbC_sp <- mtbC %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMA")
mtbC_dom_gen <- mtbC %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMA")

mtmC <- read.delim("exportdata_mtmC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtmC_dom <- mtmC %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "MMA")
mtmC_gen <- mtmC %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "MMA")
mtmC_sp <- mtmC %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "MMA")
mtmC_dom_gen <- mtmC %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "MMA")

mtsA <- read.delim("exportdata_mtsA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  filter(High.Quality == "Yes") %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtsA_dom <- mtsA %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMS/MeSH/MMPA")
mtsA_gen <- mtsA %>%
  group_by(Genus) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMS/MeSH/MMPA")
mtsA_sp <- mtsA %>%
  group_by(GenSpp) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMS/MeSH/MMPA")
mtsA_dom_gen <- mtsA %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "unclassified") %>%
  filter(Genus != "Uncultured") %>%
  filter(Genus != "uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "DMS/MeSH/MMPA")

domain <- rbind(mcrA_dom,
                cdhD_dom,
                mtrA_dom,
                mtaA_dom,
                mttC_dom,
                mtbC_dom,
                mtmC_dom,
                mtsA_dom) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("All",
                                       "Acetate",
                                       "H2/CO2",
                                       "Methanol",
                                       "TMA",
                                       "DMA",
                                       "MMA",
                                       "DMS/MeSH/MMPA"))) %>%
  mutate(Dataset = "Number of Genomes")

dom_genus <- rbind(mcrA_dom_gen,
               cdhD_dom_gen,
               mtrA_dom_gen,
               mtaA_dom_gen,
               mttC_dom_gen,
               mtbC_dom_gen,
               mtmC_dom_gen,
               mtsA_dom_gen) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("All",
                                       "Acetate",
                                       "H2/CO2",
                                       "Methanol",
                                       "TMA",
                                       "DMA",
                                       "MMA",
                                       "DMS/MeSH/MMPA"))) %>%
  mutate(Dataset = "Number of Genera")

dom_gen <- rbind(domain, dom_genus)

ggplot(dom_gen, aes(Substrate, count, fill = Domain)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, 
       y = NULL) +
  scale_fill_manual(values = c("#F8766D", "#619CFF", "grey")) +
  scale_x_discrete(limits = rev(levels(dom_gen$Substrate))) +
  coord_flip() +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))



#### All 3 methylated amines ####
# All use mtbA
mtbA <- read.delim("exportdata_mtbA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
mtbA_dom <- mtbA %>%
  group_by(Domain) %>%
  summarise(count = n())
mtbA_gen <- mtbA %>%
  group_by(Genus) %>%
  summarise(count = n())
mtbA_sp <- mtbA %>%
  group_by(GenSpp) %>%
  summarise(count = n())


#### Trimethylamine ####
TMA_mttB <- read.delim("exportdata_mttB_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
TMA_mttb_dom <- TMA_mttB %>%
  group_by(Domain) %>%
  summarise(count = n())
TMA_mttb_gen <- TMA_mttB %>%
  group_by(Genus) %>%
  summarise(count = n())
TMA_mttb_sp <- TMA_mttB %>%
  group_by(GenSpp) %>%
  summarise(count = n())

TMA_mttC <- read.delim("exportdata_mttC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
TMA_mttC_dom <- TMA_mttC %>%
  group_by(Domain) %>%
  summarise(count = n())
TMA_mttC_gen <- TMA_mttC %>%
  group_by(Genus) %>%
  summarise(count = n())
TMA_mttC_sp <- TMA_mttC %>%
  group_by(GenSpp) %>%
  summarise(count = n())

TMA_both <- TMA_mttB %>%
  filter(Genome %in% TMA_mttC$Genome) %>%
  filter(Genome %in% mtbA$Genome)
TMA_both_dom <- TMA_both %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Trimethylamine")
TMA_both_dom_gen <- TMA_both %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Trimethylamine")
TMA_both_gen <- TMA_both %>%
  group_by(Genus) %>%
  summarise(count = n())
TMA_both_sp <- TMA_both %>%
  group_by(GenSpp) %>%
  summarise(count = n())



#### Dimethylamine ####
DMA_mtbB <- read.delim("exportdata_mtbB_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
DMA_mtbB_dom <- DMA_mtbB %>%
  group_by(Domain) %>%
  summarise(count = n())
DMA_mtbB_gen <- DMA_mtbB %>%
  group_by(Genus) %>%
  summarise(count = n())
DMA_mtbB_sp <- DMA_mtbB %>%
  group_by(GenSpp) %>%
  summarise(count = n())

DMA_mtbC <- read.delim("exportdata_mtbC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
DMA_mtbC_dom <- DMA_mtbC %>%
  group_by(Domain) %>%
  summarise(count = n())
DMA_mtbC_gen <- DMA_mtbC %>%
  group_by(Genus) %>%
  summarise(count = n())
DMA_mtbC_sp <- DMA_mtbC %>%
  group_by(GenSpp) %>%
  summarise(count = n())

DMA_both <- DMA_mtbB %>%
  filter(Genome %in% DMA_mtbC$Genome) %>%
  filter(Genome %in% mtbA$Genome)
DMA_both_dom <- DMA_both %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Dimethylamine")
DMA_both_dom_gen <- DMA_both %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Dimethylamine")
DMA_both_gen <- DMA_both %>%
  group_by(Genus) %>%
  summarise(count = n())
DMA_both_sp <- DMA_both %>%
  group_by(GenSpp) %>%
  summarise(count = n())



#### Methylamine ####
MMA_mtmB <- read.delim("exportdata_mtmB_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
MMA_mtmB_dom <- MMA_mtmB %>%
  group_by(Domain) %>%
  summarise(count = n())
MMA_mtmB_gen <- MMA_mtmB %>%
  group_by(Genus) %>%
  summarise(count = n())
MMA_mtmB_sp <- MMA_mtmB %>%
  group_by(GenSpp) %>%
  summarise(count = n())

MMA_mtmC <- read.delim("exportdata_mtmC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
MMA_mtmC_dom <- MMA_mtmC %>%
  group_by(Domain) %>%
  summarise(count = n())
MMA_mtmC_gen <- MMA_mtmC %>%
  group_by(Genus) %>%
  summarise(count = n())
MMA_mtmC_sp <- MMA_mtmC %>%
  group_by(GenSpp) %>%
  summarise(count = n())

MMA_both <- MMA_mtmB %>%
  filter(Genome %in% MMA_mtmC$Genome) %>%
  filter(Genome %in% mtbA$Genome)
MMA_both_dom <- MMA_both %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methylamine")
MMA_both_dom_gen <- MMA_both %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methylamine")
MMA_both_gen <- MMA_both %>%
  group_by(Genus) %>%
  summarise(count = n())
MMA_both_sp <- MMA_both %>%
  group_by(GenSpp) %>%
  summarise(count = n())



#### Methanol ####
# mtaB, mtaC, mtaA

Met_mtaB <- read.delim("exportdata_mtaB_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
Met_mtaB_dom <- Met_mtaB %>%
  group_by(Domain) %>%
  summarise(count = n())
Met_mtaB_gen <- Met_mtaB %>%
  group_by(Genus) %>%
  summarise(count = n())
Met_mtaB_sp <- Met_mtaB %>%
  group_by(GenSpp) %>%
  summarise(count = n())

Met_mtaC <- read.delim("exportdata_mtaC_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
Met_mtaC_dom <- Met_mtaC %>%
  group_by(Domain) %>%
  summarise(count = n())
Met_mtaC_gen <- Met_mtaC %>%
  group_by(Genus) %>%
  summarise(count = n())
Met_mtaC_sp <- Met_mtaC %>%
  group_by(GenSpp) %>%
  summarise(count = n())

Met_mtaA <- read.delim("exportdata_mtaA_toFilter.txt") %>%
  mutate(Genome = Genome.Name...Sample.Name) %>%
  select(-X) %>%
  separate(Genome.Name...Sample.Name, 
           into = c("Genus", "Species", 
                    "Strain1", "Strain2", "Strain3", "Strain4", "Strain5", "Strain6", "Strain7"),
           sep = " ") %>%
  mutate(GenSpp = paste(Genus, Species), sep = "_")
Met_mtaA_dom <- Met_mtaA %>%
  group_by(Domain) %>%
  summarise(count = n())
Met_mtaA_gen <- Met_mtaA %>%
  group_by(Genus) %>%
  summarise(count = n())
Met_mtaA_sp <- Met_mtaA %>%
  group_by(GenSpp) %>%
  summarise(count = n())

Met_both <- Met_mtaB %>%
  filter(Genome %in% Met_mtaC$Genome) %>%
  filter(Genome %in% Met_mtaA$Genome)
Met_both_dom <- Met_both %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
Met_both_dom_gen <- Met_both %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
Met_both_gen <- Met_both %>%
  group_by(Genus) %>%
  summarise(count = n())
Met_both_sp <- Met_both %>%
  group_by(GenSpp) %>%
  summarise(count = n())



#### Graphs ####
domain <- rbind(TMA_both_dom,
                DMA_both_dom,
                MMA_both_dom,
                Met_both_dom) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("Methanol",
                                       "Methylamine",
                                       "Dimethylamine",
                                       "Trimethylamine"))) %>%
  mutate(Dataset = "Number of Genomes")

genus <- rbind(TMA_both_dom_gen,
                DMA_both_dom_gen,
                MMA_both_dom_gen,
                Met_both_dom_gen) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("Methanol",
                                       "Methylamine",
                                       "Dimethylamine",
                                       "Trimethylamine"))) %>%
  mutate(Dataset = "Number of Genera")

dom_gen <- rbind(domain, genus)

ggplot(dom_gen, aes(Substrate, count, fill = Domain)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, 
       y = NULL) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  coord_flip() +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))

ggplot(dom_gen, aes(Substrate, count, fill = Domain)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, 
       y = NULL) +
  scale_fill_manual(values = c("#F8766D", "#619CFF")) +
  coord_flip() +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))



#### First 2 ####
# ^ The above requires all 3 genes in each to be present. 
# What if we look at just the first 2 genes (use this)
# Revisions: split into Archaea (with mcrA), Archaea (no mcrA)
TMA_two <- TMA_mttB %>%
  filter(Genome %in% TMA_mttC$Genome) %>%
  mutate(Domain2 = Domain)
for (i in 1:nrow(TMA_two)) {
  if (TMA_two$Domain[i] == "Archaea" & TMA_two$taxon_oid[i] %in% mcrA$taxon_oid) {
    TMA_two$Domain2[i] <- "Archaea (with mcrA)"
  }
}
for (i in 1:nrow(TMA_two)) {
  if (TMA_two$Domain[i] == "Archaea" & TMA_two$taxon_oid[i] %notin% mcrA$taxon_oid) {
    TMA_two$Domain2[i] <- "Archaea (no mcrA)"
  }
}
table(TMA_two$Domain)
table(TMA_two$Domain2)
TMA_two_dom <- TMA_two %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Trimethylamine")
TMA_two_dom_gen <- TMA_two %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain2, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Trimethylamine")
TMA_two_gen <- TMA_two %>%
  group_by(Genus) %>%
  summarise(count = n())
TMA_two_sp <- TMA_two %>%
  group_by(GenSpp) %>%
  summarise(count = n())

DMA_two <- DMA_mtbB %>%
  filter(Genome %in% DMA_mtbC$Genome) %>%
  mutate(Domain2 = Domain)
for (i in 1:nrow(DMA_two)) {
  if (DMA_two$Domain[i] == "Archaea" & DMA_two$taxon_oid[i] %in% mcrA$taxon_oid) {
    DMA_two$Domain2[i] <- "Archaea (with mcrA)"
  }
}
for (i in 1:nrow(DMA_two)) {
  if (DMA_two$Domain[i] == "Archaea" & DMA_two$taxon_oid[i] %notin% mcrA$taxon_oid) {
    DMA_two$Domain2[i] <- "Archaea (no mcrA)"
  }
}
table(DMA_two$Domain)
table(DMA_two$Domain2)
DMA_two_dom <- DMA_two %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Dimethylamine")
DMA_two_dom_gen <- DMA_two %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain2, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Dimethylamine")
DMA_two_gen <- DMA_two %>%
  group_by(Genus) %>%
  summarise(count = n())
DMA_two_sp <- DMA_two %>%
  group_by(GenSpp) %>%
  summarise(count = n())

MMA_two <- MMA_mtmB %>%
  filter(Genome %in% MMA_mtmC$Genome) %>%
  mutate(Domain2 = Domain)
for (i in 1:nrow(MMA_two)) {
  if (MMA_two$Domain[i] == "Archaea" & MMA_two$taxon_oid[i] %in% mcrA$taxon_oid) {
    MMA_two$Domain2[i] <- "Archaea (with mcrA)"
  }
}
for (i in 1:nrow(MMA_two)) {
  if (MMA_two$Domain[i] == "Archaea" & MMA_two$taxon_oid[i] %notin% mcrA$taxon_oid) {
    MMA_two$Domain2[i] <- "Archaea (no mcrA)"
  }
}
table(MMA_two$Domain)
table(MMA_two$Domain2)
MMA_two_dom <- MMA_two %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methylamine")
MMA_two_dom_gen <- MMA_two %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain2, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methylamine")
MMA_two_gen <- MMA_two %>%
  group_by(Genus) %>%
  summarise(count = n())
MMA_two_sp <- MMA_two %>%
  group_by(GenSpp) %>%
  summarise(count = n())

Met_two <- Met_mtaB %>%
  filter(Genome %in% Met_mtaC$Genome) %>%
  mutate(Domain2 = Domain)
for (i in 1:nrow(Met_two)) {
  if (Met_two$Domain[i] == "Archaea" & Met_two$taxon_oid[i] %in% mcrA$taxon_oid) {
    Met_two$Domain2[i] <- "Archaea (with mcrA)"
  }
}
for (i in 1:nrow(Met_two)) {
  if (Met_two$Domain[i] == "Archaea" & Met_two$taxon_oid[i] %notin% mcrA$taxon_oid) {
    Met_two$Domain2[i] <- "Archaea (no mcrA)"
  }
}
table(Met_two$Domain)
table(Met_two$Domain2)
Met_two_dom <- Met_two %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
Met_two_dom_gen <- Met_two %>%
  filter(Genus != "Unclassified") %>%
  filter(Genus != "Uncultured") %>%
  group_by(Domain2, Genus) %>%
  summarise(count = n()) %>%
  group_by(Domain2) %>%
  summarise(count = n()) %>%
  mutate(Substrate = "Methanol")
Met_two_gen <- Met_two %>%
  group_by(Genus) %>%
  summarise(count = n())
Met_two_sp <- Met_two %>%
  group_by(GenSpp) %>%
  summarise(count = n())

domain_two <- rbind(TMA_two_dom,
                DMA_two_dom,
                MMA_two_dom,
                Met_two_dom) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("Methanol",
                                       "Methylamine",
                                       "Dimethylamine",
                                       "Trimethylamine"))) %>%
  mutate(Dataset = "Number of Genomes")

genus_two <- rbind(TMA_two_dom_gen,
               DMA_two_dom_gen,
               MMA_two_dom_gen,
               Met_two_dom_gen) %>%
  mutate(Substrate = factor(Substrate,
                            levels = c("Methanol",
                                       "Methylamine",
                                       "Dimethylamine",
                                       "Trimethylamine"))) %>%
  mutate(Dataset = "Number of Genera")

dom_gen_two <- rbind(domain_two, genus_two)
show_col(brewer_pal(palette = "Paired")(12))

pdf("Manuscript/Figure3.pdf", width = 8, height = 4)
ggplot(dom_gen_two, aes(Substrate, count, fill = Domain2)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, 
       y = NULL,
       fill = "Domain") +
  scale_fill_manual(values = c("#FB9A99", "#E31A1C", "#A6CEE3")) +
  coord_flip() +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))
dev.off()

png("Manuscript/Figure3.png", width = 8, height = 4, units = "in", res = 300)
ggplot(dom_gen_two, aes(Substrate, count, fill = Domain2)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = NULL, 
       y = NULL) +
  scale_fill_manual(values = c("#FB9A99", "#E31A1C", "#A6CEE3")) +
  coord_flip() +
  facet_wrap(~ Dataset, scales = "free_x") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12))
dev.off()

table(mtbA$Genus)



#### Methylated Amine Overlap ####
# Do same archaea use TMA, DMA and MMA
# Or are DMA, and MMA, users separate?
# Need to extract genomes with genes, and check overlap using "IMG.Genome.ID"

# TMA - mttB, mttC, mtbA, mcrA. n = 115
TMA_archaea <- TMA_mttB %>%
  filter(IMG.Genome.ID %in% TMA_mttC$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mtbA$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mcrA$IMG.Genome.ID)

# DMA - mtbB, mtbC, mtbA, mcrA. n = 113
DMA_archaea <- DMA_mtbB %>%
  filter(IMG.Genome.ID %in% DMA_mtbC$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mtbA$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mcrA$IMG.Genome.ID)

# MMA - mtmB, mtmC, mtbA, mcrA. n = 103
MMA_archaea <- MMA_mtmB %>%
  filter(IMG.Genome.ID %in% MMA_mtmC$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mtbA$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% mcrA$IMG.Genome.ID)

# Total
amines <- rbind(TMA_archaea, DMA_archaea, MMA_archaea) %>%
  group_by(IMG.Genome.ID) %>%
  slice_head(n = 1)

# Total 117
# 97 TMA/DMA/MMA
# 14 TMA/DMA
# 4 TMA/MMA
# 2 DMA/MMA
# 0 only 1

# Combos
combo_TMA_DMA_MMA <- TMA_archaea %>%
  filter(IMG.Genome.ID %in% DMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %in% MMA_archaea$IMG.Genome.ID)

combo_TMA_only <- TMA_archaea %>%
  filter(IMG.Genome.ID %notin% DMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% MMA_archaea$IMG.Genome.ID)

combo_DMA_only <- DMA_archaea %>%
  filter(IMG.Genome.ID %notin% TMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% MMA_archaea$IMG.Genome.ID)

combo_MMA_only <- MMA_archaea %>%
  filter(IMG.Genome.ID %notin% TMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% DMA_archaea$IMG.Genome.ID)

combo_TMA_DMA_only <- TMA_archaea %>%
  filter(IMG.Genome.ID %in% DMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% MMA_archaea$IMG.Genome.ID)

combo_TMA_MMA_only <- TMA_archaea %>%
  filter(IMG.Genome.ID %in% MMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% DMA_archaea$IMG.Genome.ID)

combo_DMA_MMA_only <- DMA_archaea %>%
  filter(IMG.Genome.ID %in% MMA_archaea$IMG.Genome.ID) %>%
  filter(IMG.Genome.ID %notin% TMA_archaea$IMG.Genome.ID)

