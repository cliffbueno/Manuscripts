# Figures for NWT seed microbiome/PSF manuscript
# by Cliff Bueno de Mesquita, Fall 2021
# Updated August 2022 for final revision
# Export as .png (for GoogleDoc) and .pdf (better, for submission)

#### Setup #####
# Libraries
library(plyr)
library(tidyverse)
library(V.PhyloMaker)
library(ape)
library(phytools)
library(readxl)
library(scales)

# First need plant phylogeny
sp.list <- read.csv("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/SpeciesList.csv")
tre <- phylo.maker(sp.list = sp.list, 
                   tree = GBOTB.extended, 
                   nodes = nodes.info.1, 
                   output.sp.list = TRUE, 
                   output.tree = FALSE, scenarios = "S3", 
                   r = 1)
tree <- read.newick(text = write.tree(tre$scenario.3))
plot(tree)

# Now read in data
bact <- read.csv("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/bactdata.csv")
bact$Species <- factor(bact$Species,
                       levels = c("LuzSpi", "KobMyo", "CarPyr", "TriSpi", 
                                  "DesCes", "FesBra", "GeuRos", "OxyDig",
                                  "SilAca", "EriSim"))
fung <- read.csv("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/fungdata.csv")
fung$Species <- factor(fung$Species,
                       levels = c("LuzSpi", "KobMyo", "CarPyr", "TriSpi", 
                                  "DesCes", "FesBra", "GeuRos", "OxyDig",
                                  "SilAca", "EriSim"))
names(bact) <- c("X", "ID", "Species", "Site", "Reads", "Richness", "Shannon",
                 "Axis01", "Axis02", "Dataset")
names(fung) <- c("X", "ID", "Species", "Site", "Reads", "Richness", "Shannon",
                 "Axis01", "Axis02", "Dataset")
d <- rbind(bact, fung) %>%
  mutate(Site = recode(Site,
                       "Divide" = "Navajo",
                       "Ridge" = "Saddle",
                       "Valley" = "GLV")) %>%
  mutate(Site = factor(Site, levels = c("Navajo", "GLV", "Saddle"))) %>%
  mutate(Dataset = recode(Dataset,
                          "Bacteria" = "(a) Bacteria",
                          "Fungi" = "(b) Fungi"))



#### 1 Field site (outside of R) ####
# Site map made in QGIS



#### 2 Richness 16S, ITS ####
pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/Manuscript/Figure2.pdf", 
    width = 6, height = 4)
ggplot(d, aes(x = Species, y = Richness)) +
  geom_boxplot() +
  geom_point(size = 2, aes(shape = Site)) +
  labs(x = "Species",
       y = "OTU Richness") +
  scale_x_discrete(limits = rev(levels(d$Species))) +
  facet_wrap(~ Dataset, ncol = 2) +
  coord_flip() +
  theme_bw() +
  guides("shape" = guide_legend(override.aes = list(size = 2.5))) +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.25),
        # legend.margin = margin(0.1, 0.1, 0.2, -0.3, unit = "cm"),
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, color = "black"))
dev.off()
  


#### 3 PCA 16S, ITS ####
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(d, c("Dataset", "Species"), find_hull)

text <- data.frame("Dataset" = c("(a) Bacteria", "(b) Fungi"),
                   "x" = c(0.2, -0.2),
                   "y" = c(-0.2, -0.2),
                   "label" = c("66.2%, 16.7%", "35.8%, 12.8%"))

d1 <- d %>%
  filter(Species == "LuzSpi" | Species == "KobMyo" |
           Species == "CarPyr" | Species == "DesCes" |
           Species == "GeuRos" | Species == "OxyDig" |
           Species == "SilAca" | Species == "EriSim")
d2 <- d %>% 
  filter(Species == "TriSpi" | Species == "FesBra")

pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/PCoA.pdf", width = 6, height = 3.5)
ggplot(d, aes(x = Axis01, y = Axis02, colour = Species)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Species, fill = Species), 
               alpha = 0.1, 
               show.legend = F) +
  geom_point(data = d1,
             size = 2, 
             alpha = 0.4, 
             aes(shape = Site, x = Axis01, y = Axis02)) +
  geom_point(data = d2,
             size = 3.5, 
             aes(shape = Site, x = Axis01, y = Axis02)) +
#  geom_text(data = text, 
#            aes(x = x, y = y, label = label),
#            size = 3, 
#            inherit.aes = F) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
#  scale_colour_manual(values = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
#                                 "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD")) +
  labs(x = "PC1",
       y = "PC2") +
  facet_wrap(~ Dataset, ncol = 2) +
  guides("shape" = guide_legend(override.aes = list(size = 2.5))) +
  guides("colour" = guide_legend(override.aes = list(size = c(2.5)))) +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.25),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.35, "cm"),
        # legend.margin = margin(-0.1, 0.1, 0, -0.3, unit = "cm"),
        axis.title = element_text(face = "bold", size = 14, color = "black"), 
        axis.text = element_text(size = 12, color = "black"), 
        strip.text = element_text(size = 12, color = "black"))
dev.off()



#### 4 Germination  ####
germ_s <- read_excel("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/20211023_alpinegerm.xlsx",
                   sheet = 1) %>%
  filter(treat != "slurry") %>%
  arrange(plant, treat, germ) %>%
  select(plant, treat, germ) %>%
  mutate(dataset = "Germination probability")

germ_s_sum <- germ_s %>%
  group_by(plant, treat) %>%
  summarise(n = n(),
            yes = sum(germ),
            no = n - yes,
            perc = yes/n*100)

tri <- subset(germ_s, plant == "TriSpi")
m <- glm(germ ~ 0 + treat, family = binomial, data = tri)
summary(m)
plogis(coef(m))
plogis(confint(m))

fes <- subset(germ_s, plant == "FesBra")
m1 <- glm(germ ~ 0 + treat, family = binomial, data = fes)
summary(m1)
plogis(coef(m1))
plogis(confint(m1))

ci <- read.csv("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/germ_confints.csv")

germ_t <- read_excel("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/20211023_alpinegerm.xlsx",
                     sheet = 1) %>%
  filter(treat != "slurry") %>%
  filter(germq == "Yes") %>%
  arrange(plant, treat, germday) %>%
  select(plant, treat, germday) %>%
  mutate(dataset = "Days to germination") %>%
  `colnames<-`(c("plant", "treat", "germ", "dataset"))

germ <- rbind(germ_s, germ_t)
germ$dataset <- factor(germ$dataset,
                       levels = c("Germination probability",
                                  "Days to germination"))

facet_names <- c("FesBra" = "(a) F. brachyphylla",
                 "TriSpi" = "(b) T. spicatum",
                 "Germination probability" = "Germination probability",
                 "Days to germination" = "Days to germination")

pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/Manuscript/Figure4.pdf", width = 5, height = 4)
ggplot(germ, aes(treat, germ)) +
  geom_boxplot(data = subset(germ, dataset == "Days to germination"),
               outlier.shape = NA) +
  geom_jitter(data = subset(germ, dataset == "Days to germination"),
              size = 2, alpha = 0.5, height = 0, width = 0.25) +
  geom_segment(aes(xend = treat, y = lo, yend = hi),
               data = ci, size=3, col = "blue", alpha = 0.3) + 
  geom_point(aes(x = treat, y = prob), data = ci, col = "black", size = 3) +
  labs(x = "Treatment",
       y = NULL) +
  facet_grid(dataset ~ plant, scales = "free_y", 
             labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"), 
        strip.text.x = element_text(size = 10, color = "black", face = "italic"),
        strip.text.y = element_text(size = 10, color = "black"))
dev.off()



#### 5 Biomass ####
biomass <- read_excel("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/20211023_alpinepsf.xlsx",
                      sheet = 1) %>%
  mutate_if(is.character, as.factor)
biomass$soiltype <- NA
for (i in 1:nrow(biomass)) {
  ifelse(biomass$plantsp[i] == biomass$sourcesp[i],
         biomass$soiltype[i] <- "Conspecific",
         biomass$soiltype[i] <- "Heterospecific")
}
pdf("~/Desktop/CU/2Research/Seeds/Figure5.pdf", width = 6, height = 4)
ggplot(biomass, aes(plantsp, tmass, colour = soiltype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.75, 
             position = position_jitterdodge(jitter.width = 0.2)) +
  labs(x = "Plant",
       y = "Total biomass (g)",
       colour = "Soil Source") +
  scale_colour_manual(values = c("#F8766D", "#619CFF")) +
  facet_wrap(~ plot, ncol = 2) +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.5, "cm"),
        axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"), 
        strip.text = element_text(size = 12, color = "black"))
dev.off()

# Also plot all on one to see highest
biomass <- biomass %>%
  mutate(plot2 = plot) %>%
  separate(plot2, into = c("Junk", "densrich"), sep = " ") %>%
  select(-Junk) %>%
  mutate(plotsoiltype = paste(densrich, soiltype))

fes_biomass <- subset(biomass, plantsp == "FesBra") %>%
  mutate(densrich = factor(densrich, 
                           levels = c("Low/Low",
                                      "Low/High", 
                                      "High/Low",
                                      "High/High")))

pdf("~/Desktop/CU/2Research/Seeds/S3FesPPT.pdf", width = 5, height = 4)
ggplot(fes_biomass, aes(reorder(plotsoiltype, tmass, mean), tmass, 
                    shape = soiltype, colour = densrich)) +
  geom_boxplot(aes(reorder(plotsoiltype, tmass, mean), tmass, colour = densrich),
               show.legend = F, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.1) +
  labs(x = "Soil Type",
       y = "Total biomass (g)",
       colour = "Plant density/richness",
       shape = "Soil source") +
  scale_colour_manual(values = c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080")) +
  ggtitle(label = "(a) F. brachyphylla") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.5, "cm"),
        axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black",
                                   angle = 45, hjust = 1),
        strip.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.7), units = "cm"))
dev.off()

tri_biomass <- subset(biomass, plantsp == "TriSpi") %>%
  mutate(densrich = factor(densrich, 
                           levels = c("Low/Low",
                                      "Low/High", 
                                      "High/Low",
                                      "High/High")))

pdf("~/Desktop/CU/2Research/Seeds/S3TriPPT.pdf", width = 5, height = 4)
ggplot(tri_biomass, aes(reorder(plotsoiltype, tmass, mean), tmass, 
                        shape = soiltype, colour = densrich)) +
  geom_boxplot(aes(reorder(plotsoiltype, tmass, mean), tmass, colour = densrich),
               show.legend = F, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.1) +
  labs(x = "Soil Type",
       y = "Total biomass (g)",
       colour = "Plant density/richness",
       shape = "Soil source") +
  scale_colour_manual(values = c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080")) +
  ggtitle(label = "(b) T. spicatum") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.5, "cm"),
        axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black",
                                   angle = 45, hjust = 1),
        strip.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.7), units = "cm"))
dev.off()

# Continuous
ggplot(fes_biomass, aes(richness, tmass, 
                        colour = soiltype)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Plant richness",
       y = "Total biomass (g)",
       colour = "Soil source") +
  ggtitle(label = "(a) F. brachyphylla") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.5, "cm"),
        axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), units = "cm"))
ggplot(tri_biomass, aes(richness, tmass, 
                        colour = soiltype)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm") +
  labs(x = "Plant richness",
       y = "Total biomass (g)",
       colour = "Soil source") +
  ggtitle(label = "(b) T. spicatum") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.5, "cm"),
        axis.title = element_text(face = "bold", size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        strip.text = element_text(size = 12, color = "black"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.5), units = "cm"))



#### S1 Soil Microbiome ####
# Extract from Dorota

### S2 Taxa Abundance ####
# Heatmaps made in their respective scripts and combined in PPT

#### S3 Seed mass ####
# Richness and seed weight
keep <- data.frame(Keep = c("SILACA", "OXYDIG", "ERISIM", "GEUROS", "TRISPI",
                            "KOBMYO", "FESBRA", "CARPYR", "DESCES"))
sw <- read.csv("~/Desktop/OneDrive - UCB-O365/CU/7Niwot/SeedWeights.csv") %>%
  filter(Species %in% keep$Keep) %>%
  group_by(Species) %>%
  summarise(Weight = round(mean(Weight), digits = 3)) %>%
  ungroup()

wt_rich_b <- bact %>%
  group_by(Species) %>%
  summarise(Rich = mean(Richness)) %>%
  ungroup() %>%
  filter(Species != "LuzSpi") %>%
  mutate(Sp = recode(Species,
                     "SilAca" = "SILACA",
                     "OxyDig" = "OXYDIG", 
                     "EriSim" = "ERISIM", 
                     "GeuRos" = "GEUROS", 
                     "TriSpi" = "TRISPI",
                     "KobMyo" = "KOBMYO", 
                     "FesBra" = "FESBRA",
                     "CarPyr" = "CARPYR", 
                     "DesCes" = "DESCES")) %>%
  left_join(., sw, by = c("Sp" = "Species")) %>%
  mutate(Dataset = "(a) Bacteria")

wt_rich_f <- fung %>%
  group_by(Species) %>%
  summarise(Rich = mean(Richness)) %>%
  ungroup() %>%
  filter(Species != "LuzSpi") %>%
  mutate(Sp = recode(Species,
                     "SilAca" = "SILACA",
                     "OxyDig" = "OXYDIG", 
                     "EriSim" = "ERISIM", 
                     "GeuRos" = "GEUROS", 
                     "TriSpi" = "TRISPI",
                     "KobMyo" = "KOBMYO", 
                     "FesBra" = "FESBRA",
                     "CarPyr" = "CARPYR", 
                     "DesCes" = "DESCES")) %>%
  left_join(., sw, by = c("Sp" = "Species")) %>%
  mutate(Dataset = "(b) Fungi")

wt_rich <- rbind(wt_rich_b, wt_rich_f)

m <- lm(Rich ~ Weight, data = wt_rich_b)
summary(m)
m <- lm(Rich ~ poly(Weight, 2, raw = T), data = wt_rich_b)
summary(m)
m <- lm(Rich ~ Weight, data = wt_rich_f)
summary(m)
m <- lm(Rich ~ poly(Weight, 2, raw = T), data = wt_rich_f)
summary(m)

pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Seeds/Manuscript/FigureS3.pdf", 
    width = 6, height = 3)
ggplot(wt_rich, aes(x = Weight, y = Rich, colour = Species, shape = Species)) +
  geom_point(size = 2) +
  stat_smooth(data = subset(wt_rich, Dataset == "(a) Bacteria"),
              aes(x = Weight, y = Rich),
              method = "lm", 
              formula = y ~ x + I(x^2), 
              size = 0.1, se = F,
              inherit.aes = F) +
  labs(x = "Seed mass (g/100 seeds)",
       y = "OTU richness") +
  scale_colour_manual(values = viridis_pal()(10)[2:10]) +
  scale_shape_manual(values = c(1,2,3,4,5,6,7,8,9)) +
  facet_wrap(~ Dataset, ncol = 2) +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black", size = 0.25),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        legend.key.height = unit(0.35, "cm"),
        axis.title = element_text(face = "bold", size = 14, color = "black"), 
        axis.text = element_text(size = 12, color = "black"), 
        strip.text = element_text(size = 12, color = "black"))
dev.off()


  
