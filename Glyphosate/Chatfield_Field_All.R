# Alpha and Beta Diversity Analysis
# 16S and 18S, Herbicide Project, Chatfield Field Sampling 1-4
# Cliff Bueno de Mesquita, October 2019

# Libraries
library(vegan)
library(plyr)
library(tidyverse)
library(RVAideMemoire)
library(car)
library(cowplot)
library(phyloseq)
library(picante)
library(zCompositions)
library(compositions)
library(lme4)
library(emmeans)
theme_set(theme_bw())



########################################### 16S ##########################################
# Read in Data
setwd("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide/Fall2019/16S/ChatfieldFieldAll/")
meta16 <- read.csv("ChatfieldAll_MappingFile_16S_Rare.csv")
bact_asv <- read.csv("ASV_table_rare.csv")
bact_asv$X.SampleID == meta16$X.SampleID
# Remove Outlier
meta16 <- subset(meta16, X.SampleID != "CF_III_116_23")
bact_asv <- subset(bact_asv, X.SampleID != "CF_III_116_23")
bact_asv$X.SampleID == meta16$X.SampleID
row.names(bact_asv) <- meta16$X.SampleID
bact_asv <- bact_asv[,-1]
meta16$Herbicide <- as.factor(meta16$Herbicide)
meta16$Time <- NA
for (i in 1:nrow(meta16)) {
  if (meta16$Dataset[i] == "Chatfield1") {
    meta16$Time[i] <- "Jun2018"
  } else {
    if (meta16$Dataset[i] == "Chatfield2") {
      meta16$Time[i] <- "Sep2018"
    } else {
      if (meta16$Dataset[i] == "Chatfield3") {
        meta16$Time[i] <- "Apr2019"
      } else {
        meta16$Time[i] <- "Jun2019"
      }
    }
  }
}
meta16$Time <- as.factor(meta16$Time)
meta16$Time <- factor(meta16$Time, levels = c("Jun2018","Sep2018","Apr2019","Jun2019"))
meta16 <- meta16 %>%
  separate(X.SampleID, into = c("CF", "Rest"), sep = "CF", remove = F) %>%
  separate(Rest, into = c("Plot1", "Plot2", "Subplot", "Subplot2"),
           sep = "_", remove = F)
meta16$Plot1[36:nrow(meta16)] <- meta16$Subplot[36:nrow(meta16)]
meta16$Plot2[36:nrow(meta16)] <- meta16$Subplot2[36:nrow(meta16)]
meta16$Subplot <- NULL
meta16$Subplot2 <- NULL
meta16$Rest <- NULL
meta16$CF <- NULL
meta16$Replicate <- paste(meta16$Plot1, meta16$Plot2, sep = "_")

# Contrasts are important for Type III ANOVA
# Default is treatment. Define Helmert here
my.contrasts <- list(Herbicide = "contr.Helmert", Dataset = "contr.Helmert")



### Alpha Diversity (ASV Richness, ANOVA, Boxplot)
# Calculate number of ASV in each sample
meta16$bact_Rich <- specnumber(bact_asv)
# After importing the phyloseq objects, calculate PD
pd16 <- pd(t(tot@otu_table@.Data), tot@phy_tree)
meta16$X.SampleID == rownames(pd16)
meta16$bact_PD <- pd16$PD
# Run ANOVA Model, first check homogeneity of variance
leveneTest(bact_Rich ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)
leveneTest(bact_Rich ~ Dataset, data = meta16) # Variance homogenous (p>0.05)
m1 <- aov(bact_Rich ~ Herbicide * Dataset, data = meta16)
summary(m1) # Effect of Herbicide and Dataset (time)
Anova(m1, type = "II") # Herbicide marginal, Dataset
Anova(m1, type = "III") # NSD
shapiro.test(m1$residuals) # Residuals not normally distributed but not too bad
TukeyHSD(m1) # 5 and 0
m1 <- lmer(bact_Rich ~ Herbicide * Dataset + (1|Plot1),
          data = meta16)
Anova(m1, type = "III") # NSD
m1 <- aov(bact_Rich ~ Herbicide * Dataset + Plot1, data = meta16)
summary(m1)

# Correct model...
meta16$Plot1 <- as.factor(meta16$Plot1)
meta16$Replicate <- as.factor(meta16$Replicate)
meta16$RepNum <- as.factor(meta16$RepNum)
xtabs(~ Plot1 + Replicate, meta16)

m1 <- lmer(bact_Rich ~ Herbicide*Dataset + (1|Plot1),
           data = meta16)
Anova(m1, type = "II")
Anova(m1, type = "III")

m2 <- lmer(bact_Rich ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:Replicate),
           data = meta16)
Anova(m2, type = "III")
Anova(m2, type = "II")
emmeans(m2, list(pairwise ~ Herbicide), adjust = "tukey")
summary(glht(m2,linfct=mcp(Herbicide="Tukey")))

m3 <- lmer(bact_Rich ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:RepNum),
           data = meta16)
Anova(m3, type = "III")

m4 <- lmer(bact_Rich ~ Herbicide+Dataset + (1|Plot1/Replicate),
           data = meta16)
Anova(m4, type = "II")
emmeans(m4, list(pairwise ~ Herbicide), adjust = "tukey")
summary(glht(m4,linfct=mcp(Herbicide="Tukey")))

m5 <- lmer(bact_Rich ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:Replicate),
           contrasts = my.contrasts,
           data = meta16)
Anova(m5, type = "III")

# Use Type II or Type III?
# https://mcfromnz.wordpress.com/2011/03/02/anova-type-iiiiii-ss-explained/
# Use Type II if int. not sig.
# https://towardsdatascience.com/anovas-three-types-of-estimating-sums-of-squares-don-t-make -the-wrong-choice-91107c77a27a
# Use Type III even if sig.
# https://stats.stackexchange.com/questions/60362/choice-between-type-i-type-ii-or-type-iii-anova
# 2 answers here against Type III

# Look at means by treatment
aggregate(meta16$bact_Rich, by = list(meta16$Herbicide), FUN = mean)
# Since we have two factors, let's do means and se satellites with Herbicide on x axis and four different colours of the time points
sum16 <- ddply(meta16, c("Herbicide","Time"), summarise,
               meanRich <- mean(bact_Rich),
               seRich = se(bact_Rich))
colnames(sum16) <- c("Herbicide","Time","meanRich","seRich")
pd <- position_dodge(0.3)
g1<-ggplot() +
  geom_boxplot(data = meta16,aes(Herbicide,bact_Rich), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta16,aes(Herbicide,bact_Rich, colour = Time),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum16, 
                aes(x = Herbicide, ymin=meanRich-seRich, ymax=meanRich+seRich,colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum16, aes(Herbicide, meanRich, colour = Time),
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 1000, label = "Herbicide p = 0.11\nTime p = 0.03"),
            size = 4, colour = "black") +
#  geom_text(aes(x = 1, y = 900, label = "a"),
#            size = 4, colour = "black") +
#  geom_text(aes(x = 2, y = 900, label = "ab"),
#            size = 4, colour = "black") +
#  geom_text(aes(x = 3, y = 900, label = "b"),
#            size = 4, colour = "black") +
  labs(x = NULL,
       y = "ASV Richness",
       title = "a) 16S Richness") +
  scale_colour_discrete(labels = c("Jun2018","Sep2018","Apr2019","Jun2019")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title = element_text(face="bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
g1

# Now Phylogenetic Diversity
leveneTest(bact_PD ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)
leveneTest(bact_PD ~ Dataset, data = meta16) # Variance homogenous (p>0.05)
m1pd <- aov(bact_PD ~ Herbicide * Dataset, data = meta16)
summary(m1pd) # Effect of Herbicide and Dataset (time) 
shapiro.test(m1pd$residuals) # Residuals normally distributed
TukeyHSD(m1pd) # 5 and 0
Anova(m1pd, type = "II") # Herbicide and Dataset
Anova(m1pd, type = "III") # NSD
m1pd <- lmer(bact_PD ~ Herbicide * Dataset + (1|Plot1),
           data = meta16)
Anova(m1pd, type = "III") # NSD
m1pd <- aov(bact_PD ~ Herbicide * Dataset + Plot1, data = meta16)
summary(m1pd) # Effect of Herbicide and Dataset (time)

# Correct model...
m1 <- lmer(bact_PD ~ Herbicide*Dataset + (1|Plot1),
           data = meta16)
Anova(m1, type = "III")

m2 <- lmer(bact_PD ~ Herbicide*Dataset + (1|Plot1/Replicate),
           data = meta16)
Anova(m2, type = "III")
Anova(m2, type = "II")
emmeans(m2, list(pairwise ~ Herbicide), adjust = "tukey") # p = 0.06
summary(glht(m2,linfct=mcp(Herbicide="Tukey"))) # p = 0.02

m3 <- lmer(bact_PD ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:RepNum),
           data = meta16)
Anova(m3, type = "III")

m4 <- lmer(bact_PD ~ Herbicide+Dataset + (1|Plot1/Replicate),
           data = meta16)
Anova(m4, type = "II")
emmeans(m4, list(pairwise ~ Herbicide), adjust = "tukey") # p = 0.06
summary(glht(m4,linfct=mcp(Herbicide="Tukey"))) # p = 0.02

m5 <- lmer(bact_PD ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:Replicate),
           contrasts = my.contrasts,
           data = meta16)
Anova(m5, type = "III")
emmeans(m5, list(pairwise ~ Herbicide), adjust = "tukey") # p = 0.06
summary(glht(m5,linfct=mcp(Herbicide="Tukey"))) # p = 0.02

aggregate(meta16$bact_PD, by = list(meta16$Herbicide), FUN = mean)
sum16pd <- ddply(meta16, c("Herbicide","Time"), summarise,
               meanPD <- mean(bact_PD),
               sePD = se(bact_PD))
colnames(sum16pd) <- c("Herbicide","Time","meanPD","sePD")
pd <- position_dodge(0.3)
gpd<-ggplot() +
  geom_boxplot(data = meta16,aes(Herbicide,bact_PD), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta16,aes(Herbicide,bact_PD, colour = Time),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum16pd, aes(x = Herbicide, ymin=meanPD-sePD, ymax=meanPD+sePD,
                                  colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum16pd, aes(Herbicide, meanPD, colour = Time),
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 52, label = "Herbicide p = 0.02\nTime p = 0.006"),
            size = 4, colour = "black") +
  geom_text(aes(x = 1, y = 47, label = "a"),
            size = 4, colour = "black") +
  geom_text(aes(x = 2, y = 47, label = "ab"),
            size = 4, colour = "black") +
  geom_text(aes(x = 3, y = 47, label = "b"),
            size = 4, colour = "black") +
  labs(x = NULL,
       y = "Phylogenetic Diversity",
       title = "c) 16S PD") +
  scale_colour_discrete(labels = c("Jun2018","Sep2018","Apr2019","Jun2019")) +
  ylim(25,54) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title = element_text(face="bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
gpd



### Beta Diversity (Weighted Unifrac, PERMANOVA, PCoA)
biom_file <- paste("ASV_table_rare_pc.biom", sep = "")
map_file <- paste("ChatfieldAll_MappingFile_16S_Rare.txt", sep = "")
biom_otu_tax <- import_biom(biom_file, "rep_phylo.tre")
bmsd <- import_qiime_sample_data(map_file)
tot <- merge_phyloseq(biom_otu_tax, bmsd)
# PERMANOVA and PERMDISP on weighted UniFrac distance
varespec.uni <- UniFrac(tot, weighted = TRUE)
set.seed(100)
m <- adonis(varespec.uni ~ meta16$Herbicide * meta16$Dataset, permutations = 999)
m # Significant
# Herbicide F = 2.43, R2 = 0.03, p = 0.004
# Dataset F = 17.91, R2 = 0.28, p = 0.001
set.seed(100)
m <- adonis(varespec.uni ~ meta16$Herbicide * meta16$Dataset, 
            strata = meta16$Plot1, permutations = 999)
m # Significant with strata too
# Remove interaction
set.seed(100)
m <- adonis(varespec.uni ~ meta16$Herbicide + meta16$Dataset, 
            strata = meta16$Plot1, permutations = 999)
m
# Switch order
set.seed(100)
m <- adonis(varespec.uni ~ meta16$Dataset + meta16$Herbicide, 
            strata = meta16$Plot1, permutations = 999)
m

# adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(meta16, Plot1:Replicate)
set.seed(100)
m <- adonis2(varespec.uni ~ Herbicide + Dataset,
             data = meta16,
             permutations = perm,
             by = "margin")
m
set.seed(100)
m <- adonis2(varespec.uni ~ Herbicide * Dataset,
             data = meta16,
             permutations = perm,
             by = "margin")
m

pairwise.perm.manova(varespec.uni, fact = meta16$Herbicide, nperm = 999) # 5-0
pairwise.perm.manova(varespec.uni, fact = meta16$Dataset, nperm = 999) # All
m <- betadisper(varespec.uni, meta16$Herbicide)
anova(m) # Not significant p = 0.50, F = 0.71
m <- betadisper(varespec.uni, meta16$Dataset)
anova(m) # Significant p = 0.0003, F = 6.7

# Compute Principle Coordinates Analysis ordination on weighted UniFrac distance
ordu <- ordinate(tot, "PCoA", "unifrac", weighted = TRUE)
# Store the PCoA results
df <- as.data.frame(as.matrix(ordu$vectors))
df$sample <- row.names(df)
meta16$Axis.1 <- df$Axis.1
meta16$Axis.2 <- df$Axis.2
# Calculate percentage variation explained
eig2 <- ordu$values$Eigenvalues
eig2 / sum(eig2) # Look at 1st two numbers for % variation explained
# Calculate vectors for the graph for the top contributing species (from SIMPER analysis)
# In this case, all pairwise comparisons are significant. So graph the top species contributing to difference between the control and each treatment
rel16 <- decostand(bact_asv,"total")
sim16 <- simper(rel16, meta16$Herbicide)
s16 <- summary(sim16)
head(s16$`0_5`, n = 5) 
# 1, 8, 2, 3, 4
# Check just last timepoint
cf4meta <- subset(meta16, Dataset == "Chatfield4")
cf4rel <- rel16[106:141,]
rownames(cf4rel) == cf4meta$X.SampleID
sim16 <- simper(cf4rel, cf4meta$Herbicide)
s16 <- summary(sim16)
head(s16$`0_5`, n = 5) 
# Top 5 not contributing much. Check dataset
simd16 <- simper(rel16, meta16$Dataset)
sd16 <- summary(simd16)
head(sd16$Chatfield1_Chatfield2, n = 5)
head(sd16$Chatfield1_Chatfield3, n = 5)
head(sd16$Chatfield1_Chatfield4, n = 5)
# Look at 1, 8, 2, 4, 3
tot_transposed <- t(tot@otu_table)
w <- wascores(x = ordu$vectors, w = tot_transposed)
wdf <- as.data.frame(w)
wdf$Species <- rownames(w)
sub16 <- subset(wdf, Species == "ESV_1" | Species == "ESV_8" | Species == "ESV_2" | Species == "ESV_4" | Species == "ESV_3")
sub16 # Order is 8, 4, 3, 1, 2
# Careful, double check that names match ESVs!
# ESV_1 = Pyrinomonadaceae, ESV_8 = Tychonema, ESV_2 = Udaeobacter, ESV_4 = Nitrososphaeraceae, ESV_3 = Nitrososphaera
sub16$shortnames <- c("Tychonema","Nitrososphaeraceae","Nitrososphaera","Pyrinomonadaceae","Udaeobacter")
# Graph (PCoA with convex hull polygons, samples as points, top taxa as vectors and text). This was kind of messy. New idea - plot centroids of herbicide treatments over time. Removed CF_III_116_23 (big outlier, need to look into that...)
find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2),]
micro.hulls.16 <- ddply(meta16, "TimeHerb", find_hull)
ggplot(meta16, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.16,aes(colour=TimeHerb,fill=TimeHerb),
               alpha = 0.1, size = 0.25, show.legend = F) +
  geom_point(size = 3, aes(colour = TimeHerb)) +
  labs(x = "PC1: 33.45%", y = "PC2: 12.39%", title = "a) 16S Beta Diversity") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(face="bold", size = 16, vjust = 2), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold",size=16),
        plot.margin = unit(c(0,0.1,0,0.1),"cm"))
# Centroids
m <- betadisper(varespec.uni, meta16$TimeHerb)
centroids16 <- as.data.frame(m$centroids[,1:2])
centroids16$Herbicide <- as.factor(c("0","2","5","0","2","5","0","2","5","0","2","5"))
centroids16$Time <- as.factor(c("Jun2018","Jun2018","Jun2018","Sep2018","Sep2018","Sep2018",
                              "Apr2019","Apr2019","Apr2019","Jun2019","Jun2019","Jun2019"))
centroids16$Time <- factor(centroids16$Time, levels = c("Jun2018","Sep2018","Apr2019","Jun2019"))
micro.hulls.16 <- ddply(meta16, "Time", find_hull)
g2<-ggplot() +
  geom_point(data = meta16,
             aes(x = Axis.1, y = Axis.2, colour = Time, shape = Herbicide),
             size = 3, alpha = 0.2) +
  geom_point(data = centroids16,
             aes(x = PCoA1, y = PCoA2, colour = Time, shape = Herbicide),
             size = 5) +
  geom_text(aes(x = 0.11, y = 0.08, label = "Herbicide p = 0.001\nTime p = 0.001"), size = 4, colour = "black") +
  labs(x = "PC1: 33.45%", y = "PC2: 12.39%", title = "a) 16S Beta Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.x = element_text(face="bold", size = 14, vjust = 2), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(face="bold",size=14),
        plot.margin = unit(c(0,0.1,0,0.1),"cm"))
g2
#   geom_polygon(data = micro.hulls.16,
# aes(x = Axis.1, y = Axis.2, group = Timepoint),
# alpha = 0.1, size = 0.25, colour = "gray", fill = NA,show.legend = F) +

# Top 5 drivers abundance
# ESV_1 = Pyrinomonadaceae, ESV_8 = Tychonema, ESV_2 = Udaeobacter, ESV_4 = Nitrososphaeraceae, ESV_3 = Nitrososphaera
meta16$Pyrinomonadaceae <- bact_asv$ESV_1/8889
meta16$Tychonema <- bact_asv$ESV_8/8889
meta16$Udaeobacter <- bact_asv$ESV_2/8889
meta16$Nitrososphaeraceae <- bact_asv$ESV_4/8889
meta16$Nitrososphaera <- bact_asv$ESV_3/8889
leveneTest(Pyrinomonadaceae ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)
leveneTest(Tychonema ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)
leveneTest(Udaeobacter ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)
leveneTest(Nitrososphaeraceae ~ Herbicide, data = meta16) # Variance homogenous (p>0.01)
leveneTest(Nitrososphaera ~ Herbicide, data = meta16) # Variance homogenous (p>0.05)

meta16last <- subset(meta16, Dataset == "Chatfield4")
m <- aov(Pyrinomonadaceae ~ Herbicide, data = meta16last)
summary(m) # Dataset significant (F = 3.28, p = 0.02), Herbicide not
shapiro.test(mm$residuals) # Residuals not normally distributed (p<0.05)
# Run Tukey Posthoc test
TukeyHSD(mm) # All sig but 3 diff from 1 and 4
sum18 <- ddply(meta18, c("Herbicide","Time"), summarise,
               meanMort <- mean(Mortierella),
               seMort = se(Mortierella))
colnames(sum18) <- c("Herbicide","Time","meanMort","seMort")
pd <- position_dodge(0.3)
ggplot() +
  geom_boxplot(data = meta18,aes(Herbicide,Mortierella), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta18,aes(Herbicide,Mortierella, colour = Time, shape = Herbicide),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum18, aes(x = Herbicide, ymin=meanMort-seMort, ymax=meanMort+seMort,
                                  colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum18, aes(Herbicide, meanMort, colour = Time, shape = Herbicide), size = 4, position = pd) +
  geom_text(aes(x = 3, y = 0.2, label = "Herbicide p = 0.98\nTime p = 0.02"), size = 5, colour = "black") +
  labs(x = "Herbicide Applications", y = "Mortierella Relative Abundance") +
  scale_colour_discrete(labels = c("Jun2018","Sep2018","Apr2019","Jun2019")) +
  guides(shape = FALSE) +
  theme(legend.position = c(0.51, 0.925),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, vjust = -0.5),
        axis.title = element_text(face="bold", size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

# Aitchison
bac_asv_czm <- cmultRepl(bact_asv, label = 0, method = "CZM")
bac_asv_ilr <- ilr(bac_asv_czm)
bac_ailr <- as.matrix(vegdist(sd_bacs_ilr, method = "euclidean"))
diag(bac_ailr) <- NA
bac_ailr[lower.tri(bac_ailr)] <- NA
sd_euk_czm <- cmultRepl(t(euk_input_filt_rare11$data_loaded),  label = 0, method = "CZM")
sd_euk_ilr <- ilr(sd_euk_czm)
euk_ailr <- as.matrix(vegdist(sd_euk_ilr, method = "euclidean"))
diag(euk_ailr) <- NA
euk_ailr[lower.tri(euk_ailr)] <- NA



########################################## 18S ##############################################
# Read in Data
setwd("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide/Fall2019/18S/ChatfieldAll/")
meta18 <- read.csv("ChatfieldAll_MappingFile_18S_rare.csv")
euk_asv <- read.csv("ASV_table_rare.csv")
# Remove CF_III_116_23 (only 4 ASVs!)
meta18 <- subset(meta18, SampleID != "CF_III_116_23")
euk_asv <- subset(euk_asv, X.SampleID != "CF_III_116_23")
euk_asv$X.SampleID == meta18$SampleID
row.names(euk_asv) <- meta18$SampleID
euk_asv <- euk_asv[,-1]
meta18$Herbicide <- as.factor(meta18$Herbicide)
meta18$Time <- NA
for (i in 1:nrow(meta18)) {
  if (meta18$Dataset[i] == "Chatfield1") {
    meta18$Time[i] <- "Jun2018"
  } else {
    if (meta18$Dataset[i] == "Chatfield2") {
      meta18$Time[i] <- "Sep2018"
    } else {
      if (meta18$Dataset[i] == "Chatfield3") {
        meta18$Time[i] <- "Apr2019"
      } else {
        meta18$Time[i] <- "Jun2019"
      }
    }
  }
}
meta18$Time <- as.factor(meta18$Time)
meta18$Time <- factor(meta18$Time, levels = c("Jun2018","Sep2018","Apr2019","Jun2019"))
meta18 <- meta18 %>%
  separate(SampleID, into = c("CF", "Rest"), sep = "CF", remove = F) %>%
  separate(Rest, into = c("Plot1", "Plot2", "Subplot", "Subplot2"),
           sep = "_", remove = F)
meta18$Plot1[37:nrow(meta18)] <- meta18$Subplot[37:nrow(meta18)]
meta18$Plot2[37:nrow(meta18)] <- meta18$Subplot2[37:nrow(meta18)]
meta18$Subplot <- NULL
meta18$Subplot2 <- NULL
meta18$Rest <- NULL
meta18$CF <- NULL
meta18$Replicate <- paste(meta18$Plot1, meta18$Plot2, sep = "_")



#### Alpha Diversity (ASV Richness, ANOVA, Boxplot)
# Calculate number of ASV in each sample
meta18$euk_Rich <- specnumber(euk_asv)
# After importing the phyloseq objects, calculate PD
pd18 <- pd(t(tot@otu_table@.Data), tot@phy_tree)
meta18$SampleID == rownames(pd18)
meta18$euk_PD <- pd18$PD
# Run ANOVA Model, first check homogeneity of variance
leveneTest(euk_Rich ~ Herbicide, data = meta18) # Variance homogenous (p>0.05)
leveneTest(euk_Rich ~ Dataset, data = meta18) # Variance homogenous (p>0.05)
m2 <- aov(euk_Rich ~ Herbicide*Dataset, data = meta18)
summary(m2) # Dataset significant (F = 65.68, p < 0.001), Herbicide not
shapiro.test(m2$residuals) # Residuals normally distributed (p>0.05)
# Run Tukey Posthoc test
TukeyHSD(m2) # All sig but 4-3
Anova(m2, type = "II") # Dataset
Anova(m2, type = "III") # Dataset
m2 <- lmer(euk_Rich ~ Herbicide * Dataset + (1|Plot1/Replicate),
           data = meta18)
Anova(m2, type = "III") # Dataset
Anova(m2, type = "II") # Dataset
m2 <- aov(euk_Rich ~ Herbicide*Dataset + Plot1, data = meta18)
summary(m2) # Dataset significant (F = 78.51, p < 0.001), Herbicide not

# Correct model...
meta18$Plot1 <- as.factor(meta18$Plot1)
meta18$Replicate <- as.factor(meta18$Replicate)

m4 <- lmer(euk_Rich ~ Herbicide + Dataset + (1|Plot1/Replicate),
           data = meta18)
Anova(m4, type = "II")

m5 <- lmer(euk_Rich ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:Replicate),
           contrasts = my.contrasts,
           data = meta18)
Anova(m5, type = "III")

# Look at means by treatment
aggregate(meta18$euk_Rich, by = list(meta18$Dataset), FUN = mean)

# Since we have two factors, let's do means and se satellites with Herbicide on x axis and four different colours of the time points
sum18 <- ddply(meta18, c("Herbicide","Time"), summarise,
               meanRich <- mean(euk_Rich),
               seRich = se(euk_Rich))
colnames(sum18) <- c("Herbicide","Time","meanRich","seRich")
pd <- position_dodge(0.3)
g3 <- ggplot() +
  geom_boxplot(data = meta18,aes(Herbicide,euk_Rich), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta18,aes(Herbicide,euk_Rich, colour = Time),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum18,
                aes(x = Herbicide, ymin=meanRich-seRich, ymax=meanRich+seRich,colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum18, aes(Herbicide, meanRich, colour = Time),
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 385, label = "Herbicide p = 0.67\nTime p < 0.001"),
            size = 4, colour = "black") +
  labs(x = NULL,
       y = "ASV Richness",
       title = "b) 18S Richness") +
  scale_colour_discrete(labels = c("June2018","August2018","June2019","August2019")) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title = element_text(face="bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
g3

# Phylogenetic Diversity
leveneTest(euk_PD ~ Herbicide, data = meta18) # Variance homogenous (p>0.05)
leveneTest(euk_PD ~ Dataset, data = meta18) # Variance homogenous (p>0.05)
m2pd <- aov(euk_PD ~ Herbicide*Dataset, data = meta18)
summary(m2pd) # Dataset significant (F = 60.21, p < 0.001), Herbicide not
shapiro.test(m2pd$residuals) # Residuals normally distributed (p>0.05)
TukeyHSD(m2pd) # All sig but 3-1
Anova(m2pd, type = "II") # Dataset
Anova(m2pd, type = "III") # Dataset 
m2pd <- lmer(euk_PD ~ Herbicide * Dataset + (1|Plot1/Replicate),
           data = meta18)
Anova(m2pd, type = "III") # Dataset
Anova(m2pd, type = "II") # Dataset
m2pd <- aov(euk_PD ~ Herbicide*Dataset + Plot1, data = meta18)
summary(m2pd) # Dataset significant (F = 65.68, p < 0.001), Herbicide not

# Correct model
m4 <- lmer(euk_PD ~ Herbicide + Dataset + (1|Plot1/Replicate),
           data = meta18)
Anova(m4, type = "II")

m5 <- lmer(euk_PD ~ Herbicide*Dataset + (1|Plot1) + (1|Plot1:Replicate),
           contrasts = my.contrasts,
           data = meta18)
Anova(m5, type = "III")

aggregate(meta18$euk_PD, by = list(meta18$Dataset), FUN = mean)
sum18pd <- ddply(meta18, c("Herbicide","Time"), summarise,
               meanPD <- mean(euk_PD),
               sePD = se(euk_PD))
colnames(sum18pd) <- c("Herbicide","Time","meanPD","sePD")
pd <- position_dodge(0.3)
gpd2 <- ggplot() +
  geom_boxplot(data = meta18,aes(Herbicide,euk_PD), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta18,aes(Herbicide,euk_PD, colour = Time),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum18pd,
                aes(x = Herbicide, ymin=meanPD-sePD, ymax=meanPD+sePD,colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum18pd, aes(Herbicide, meanPD, colour = Time),
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 92, label = "Herbicide p = 0.49\nTime p < 0.001"),
            size = 4, colour = "black") +
  labs(x = NULL,
       y = "Phylogenetic Diversity",
       title = "d) 18S PD") +
  scale_colour_discrete(labels = c("June2018","August2018","June2019","August2019")) +
  ylim(30,95) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title = element_text(face="bold", size = 14), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
gpd2

# Alpha diversity multipanel
# 16S, 18S, ASV richness, PD
multiplot <- plot_grid(g1,g3,gpd,gpd2, ncol=2, nrow=2, align = "hv")
multiplot 
# 10.52 x 5.84

# Mortierella abundance
# Add up all Mortierella ASVs, then test for effects of time and herbicide
meta18$Mortierella<-euk_asv$ESV_1+euk_asv$ESV_1146+euk_asv$ESV_1151+euk_asv$ESV_12+euk_asv$ESV_1374+euk_asv$ESV_1584+euk_asv$ESV_1594+euk_asv$ESV_1942+euk_asv$ESV_1961+euk_asv$ESV_1995+euk_asv$ESV_213+euk_asv$ESV_255+euk_asv$ESV_445+euk_asv$ESV_462+euk_asv$ESV_771+euk_asv$ESV_794+euk_asv$ESV_806+euk_asv$ESV_847+euk_asv$ESV_851+euk_asv$ESV_874
meta18$Mortierella <- meta18$Mortierella/5226
leveneTest(Mortierella ~ Herbicide, data = meta18) # Variance homogenous (p>0.05)
leveneTest(Mortierella ~ Dataset, data = meta18) # Variance homogenous (p>0.05)
mm <- aov(Mortierella ~ Herbicide*Dataset, data = meta18)
summary(mm) # Dataset significant (F = 3.28, p = 0.02), Herbicide not
shapiro.test(mm$residuals) # Residuals not normally distributed (p<0.05)
# Run Tukey Posthoc test
TukeyHSD(mm) # All sig but 3 diff from 1 and 4
sum18 <- ddply(meta18, c("Herbicide","Time"), summarise,
               meanMort <- mean(Mortierella),
               seMort = se(Mortierella))
colnames(sum18) <- c("Herbicide","Time","meanMort","seMort")
pd <- position_dodge(0.3)
ggplot() +
  geom_boxplot(data = meta18,aes(Herbicide,Mortierella), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = meta18,aes(Herbicide,Mortierella, colour = Time, shape = Herbicide),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sum18, aes(x = Herbicide, ymin=meanMort-seMort, ymax=meanMort+seMort,
                                  colour = Time),
                width=.4, position=pd, size =0.9) +
  geom_point(data = sum18, aes(Herbicide, meanMort, colour = Time, shape = Herbicide), size = 4, position = pd) +
  geom_text(aes(x = 3, y = 0.2, label = "Herbicide p = 0.98\nTime p = 0.02"), size = 5, colour = "black") +
  labs(x = "Herbicide Applications", y = "Mortierella Relative Abundance") +
  scale_colour_discrete(labels = c("Jun2018","Sep2018","Apr2019","Jun2019")) +
  guides(shape = FALSE) +
  theme(legend.position = c(0.51, 0.925),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black"),
        legend.key = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16, vjust = -0.5),
        axis.title = element_text(face="bold", size = 16), 
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))



### Beta Diversity (Weighted Unifrac, PERMANOVA, PCoA)
# Import files with phyloseq package
biom_file <- paste("ASV_table_rare.biom", sep = "")
map_file <- paste("ChatfieldAll_MappingFile_18S_rare.txt", sep = "")
biom_otu_tax <- import_biom(biom_file, "rep_phylo.tre")
bmsd <- import_qiime_sample_data(map_file)
tot <- merge_phyloseq(biom_otu_tax, bmsd)
tot <- subset_samples(tot, sample_names(tot) != "CF_III_116_23")
# PERMANOVA and PERMDISP on weighted UniFrac distance
varespec.uni <- UniFrac(tot, weighted = TRUE)
set.seed(100)
m <- adonis(varespec.uni ~ meta18$Herbicide * meta18$Dataset, permutations = 999)
m
set.seed(100)
m <- adonis(varespec.uni ~ meta18$Herbicide * meta18$Dataset, 
            strata = meta18$Plot1, permutations = 999)
m
# Herbicide p = 0.14, F = 1.28, R2 = 0.02
# Time p = 0.001, F = 7.33, R2 = 0.14

# adonis2
perm <- how(nperm = 999)
setBlocks(perm) <- with(meta18, Plot1:Replicate)
set.seed(100)
m <- adonis2(varespec.uni ~ Herbicide + Dataset,
             data = meta18,
             permutations = perm,
             by = "margin")
m

set.seed(100)
m <- adonis2(varespec.uni ~ Herbicide * Dataset,
             data = meta18,
             permutations = perm,
             by = "margin")
m

pairwise.perm.manova(varespec.uni, fact = meta18$Herbicide, nperm = 999) # NSD
pairwise.perm.manova(varespec.uni, fact = meta18$Dataset, nperm = 999) # All significant

m <- betadisper(varespec.uni, meta18$Herbicide)
anova(m) # Not significant p = 0.48, F = 0.73
m <- betadisper(varespec.uni, meta18$Dataset)
anova(m) # Significant p = 0.008, F = 4.1

# Compute Principle Coordinates Analysis ordination on weighted UniFrac distance
ordu <- ordinate(tot, "PCoA", "unifrac", weighted = TRUE)
# Store the PCoA results
df <- as.data.frame(as.matrix(ordu$vectors))
df$sample <- row.names(df)
meta18$Axis.1 <- df$Axis.1
meta18$Axis.2 <- df$Axis.2
# Calculate percentage variation explained
eig2 <- ordu$values$Eigenvalues
eig2 / sum(eig2) # Look at 1st two numbers for % variation explained
# Calculate hulls for the graph
find_hull <- function(df) df[chull(df$Axis.1, df$Axis.2),]
micro.hulls.18 <- ddply(meta18, "Treatment", find_hull)
# Calculate vectors for the graph for the top contributing species (from SIMPER analysis)
# In this case, all pairwise comparisons are significant. So graph the top species contributing to difference between the control and each treatment
rel18 <- decostand(euk_asv,"total")
sim18 <- simper(rel18, meta18$Dataset)
s18 <- summary(sim18)
head(s18$Chatfield1_Chatfield2, n = 5) 
head(s18$Chatfield1_Chatfield3, n = 5) 
head(s18$Chatfield1_Chatfield4, n = 5)
# 2, 1, 6, 4 consistent contributers.
tot_transposed <- t(tot@otu_table)
w <- wascores(x = ordu$vectors, w = tot_transposed)
wdf <- as.data.frame(w)
wdf$Species <- rownames(w)
sub18 <- subset(wdf,Species=="ESV_2"|Species =="ESV_1"|Species =="ESV_6"|Species=="ESV_4")
sub18
# Order is 6, 4, 2, 1
# ESV_6 = Hyalorbilia , ESV_4 = Sclerotinia , ESV_2 = Rhizophlyctis, ESV_1 = Mortierella
sub18$shortnames <- c("Hyalorbilia","Sclerotinia","Rhizophlyctis","Mortierella")
# Graph (PCoA with convex hull polygons, samples as points, top taxa as vectors and text)
ggplot(meta18, aes(Axis.1, Axis.2)) +
  geom_polygon(data = micro.hulls.18,aes(colour=Treatment,fill=Treatment),
               alpha = 0.1, size = 0.25, show.legend = F) +
  geom_point(size = 3, aes(colour = Treatment)) +
  geom_segment(data=sub18,aes(x=0,xend=Axis.1,y=0,yend=Axis.2),
               arrow = arrow(length = unit(0.5, "cm")),colour="grey",
               inherit.aes = FALSE) + 
  geom_text(data=sub18,aes(x=Axis.1,y=Axis.2,label=shortnames),
            size=4,colour = "blue") +
  labs(x = "PC1: 17.54%", y = "PC2: 13.89%", title = "d) 18S Beta Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(face="bold", size = 16, vjust = 5), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold",size=16),
        plot.margin = unit(c(0,0.1,0,0.1),"cm"))
# Too messy, focus on centroids
# Centroids
m <- betadisper(varespec.uni, meta18$TimeHerb)
centroids18 <- as.data.frame(m$centroids[,1:2])
centroids18$Herbicide <- as.factor(c("0","2","5","0","2","5","0","2","5","0","2","5"))
centroids18$Time <- as.factor(c("Jun2018","Jun2018","Jun2018","Sep2018","Sep2018","Sep2018",
                              "Apr2019","Apr2019","Apr2019","Jun2019","Jun2019","Jun2019"))
centroids18$Time <- factor(centroids18$Time, levels = c("Jun2018","Sep2018","Apr2019","Jun2019"))
micro.hulls.18 <- ddply(meta18, "Time", find_hull)
g4<-ggplot() +
  geom_point(data = meta18,
             aes(x = Axis.1, y = Axis.2, colour = Time, shape = Herbicide),
             size = 3, alpha = 0.2) +
  geom_point(data = centroids18,
             aes(x = PCoA1, y = PCoA2, colour = Time, shape = Herbicide),
             size = 5) +
  geom_text(aes(x = 0.13, y = 0.25, label = "Herbicide p = 0.001\nTime p = 0.001"), 
            size = 4, colour = "black") +
  labs(x = "PC1: 17.54%", y = "PC2: 13.89%", title = "b) 18S Beta Diversity") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.x = element_text(face="bold", size = 14, vjust = 2), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(face="bold",size=14),
        plot.margin = unit(c(0,0.1,0,0.1),"cm"))
g4



#################################### Multipanel Graphs ######################################
# Alpha Diversity
figure3 <- plot_grid(g1,g3,gpd,gpd2, ncol=2, nrow=2, align = "hv")
figure3
# 10.63 x 5.94
pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide/Figure3_forPPT.pdf",
    width = 9, height = 6)
figure3
dev.off()

# Just beta
figure4 <- plot_grid(g2,g4, ncol=1, nrow=2, align = "hv")
figure4
# 6.01 x 6.32
pdf("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide/Figure4_forPPT.pdf",
    width = 6, height = 6)
figure4
dev.off()
