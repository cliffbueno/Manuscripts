# Analysis of soil bacterial communities and enzyme activity in moss versus vascular plant dominated plots in the Green Lakes Valley (Andrew King) and Soddie (Joey Knelman) at the Niwot Ridge LTER site
# By Cliff Bueno de Mesquita, Fall 2016
# Now published in Soil Biology and Biochemistry
# Did analyses on multiple taxonomic levels, but wrote paper for order level

# Setup
library(permute)
library(lattice)
library(vegan)
library(plyr)
library(ggplot2)
library(lmodel2)
source("~/Desktop/Functions/Multiplot.R")
setwd("~/Desktop/CU/2Research/Moss")
sub <- read.csv("2007data.csv", header = TRUE)
meta <- read.csv("Plots.csv", header = TRUE)
bactP <- t(read.table("moss_tax_L2.txt", sep="\t", header=T, row.names=1))
bactC <- t(read.table("moss_tax_L3.txt", sep="\t", header=T, row.names=1))
bactO <- t(read.table("moss_tax_L4.txt", sep="\t", header=T, row.names=1))
bactF <- t(read.table("moss_tax_L5.txt", sep="\t", header=T, row.names=1))
bactG <- t(read.table("moss_tax_L6.txt", sep="\t", header=T, row.names=1))
phylum <- as.data.frame(bactP[1:20,])
class <- as.data.frame(bactC[1:20,])
order <- as.data.frame(bactO[1:20,])
family <- bactF[1:20,]
genus <- bactG[1:20,]
# Add ratios to the dataset 
sub$bgnag <- sub$BG/sub$NAG
sub$nagphos <- sub$NAG/sub$PHOS
sub$bgphos <- sub$BG/sub$PHOS

############################## 2007 Microbial Community ####################################
# Phylum Level
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(phylum, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0.13
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
sub$Axis01 <- varespec.nmds$points[,1]
sub$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(sub, "Class", find_hull)
ggplot(sub, aes(Axis01, Axis02), colour=Class) +
    geom_polygon(data = micro.hulls, aes(colour=Class, fill=Class), 
        alpha = 0.1) +
    geom_point(aes(fill=Class), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ sub$Class, permutations=999)
type.permanova #Significant
#Taxa driving differences
phylum.pvalues<-apply(phylum, 2, function(x){wilcox.test(x~sub$Class)$p.value})
phylum.pvalues.fdr<-p.adjust(phylum.pvalues, method="fdr")
phylumPs<-data.frame(phylum.pvalues)
sigphyla<-subset(phylumPs,phylum.pvalues<0.05)
sigphyla

# Class Level
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(class, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0.15
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
sub$Axis01 <- varespec.nmds$points[,1]
sub$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(sub, "Class", find_hull)
ggplot(sub, aes(Axis01, Axis02), colour=Class) +
    geom_polygon(data = micro.hulls, aes(colour=Class, fill=Class), 
                 alpha = 0.1) +
    geom_point(aes(fill=Class), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ sub$Class, permutations=999)
type.permanova #Significant
#Taxa driving differences
class.pvalues<-apply(class, 2, function(x){wilcox.test(x~sub$Class)$p.value})
class.pvalues.fdr<-p.adjust(class.pvalues, method="fdr")
classPs<-data.frame(class.pvalues)
sigclasses<-subset(classPs,class.pvalues<0.05)
sigclasses

# Order level
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(order, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0.14
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
sub$Axis01 <- varespec.nmds$points[,1]
sub$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(sub, "Class", find_hull)
ggplot(sub, aes(Axis01, Axis02), colour=Class) +
    geom_polygon(data = micro.hulls, aes(colour=Class, fill=Class), 
                 alpha = 0.1) +
    geom_point(aes(fill=Class), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ sub$Class, permutations=999)
type.permanova #Significant
#Taxa driving differences
order.pvalues<-apply(order, 2, function(x){wilcox.test(x~sub$Class)$p.value})
order.pvalues.fdr<-p.adjust(order.pvalues, method="fdr")
orderPs<-data.frame(order.pvalues)
sigorders<-subset(orderPs,order.pvalues<0.05)
sigorders
bray <- as.matrix(varespec.bray)
# write.csv(bray, "Bray_Curtis_Order.csv")

# Specific Comparisons
wilcox.test(class$`k__Bacteria;p__Acidobacteria;c__[Chloracidobacteria]` ~ sub$Class)
wilcox.test(order$`k__Bacteria;p__Acidobacteria;c__[Chloracidobacteria];o__RB41` ~ sub$Class)
wilcox.test(phylum$`k__Bacteria;p__Verrucomicrobia` ~ sub$Class)
wilcox.test(class$`k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria]` ~ sub$Class)
wilcox.test(order$`k__Bacteria;p__Verrucomicrobia;c__[Spartobacteria];o__[Chthoniobacterales]` ~ sub$Class)

######################## Time Series Microbial Community Analysis ##########################
Tphylum <- bactP[c(8,13,14,17,21:24),]
Tclass <- bactC[c(8,13,14,17,21:24),]
Torder <- bactO[c(8,13,14,17,21:24),]
Tmeta <- meta[c(1:8),]
Tmeta$Year <- as.factor(Tmeta$Year)

# Phylum Level
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(Tphylum, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
Tmeta$Axis01 <- varespec.nmds$points[,1]
Tmeta$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(Tmeta, "Year", find_hull)
ggplot(Tmeta, aes(Axis01, Axis02), colour=Year) +
    geom_polygon(data = micro.hulls, aes(colour=Year, fill=Year), 
                 alpha = 0.1) +
    geom_point(aes(fill=Year), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ Tmeta$Year, permutations=999)
type.permanova #Significant
#Taxa driving differences
Tphylum.pvalues<-apply(Tphylum, 2, function(x){wilcox.test(x~Tmeta$Year)$p.value})
Tphylum.pvalues.fdr<-p.adjust(Tphylum.pvalues, method="fdr")
TphylumPs<-data.frame(Tphylum.pvalues)
sigTphyla<-subset(TphylumPs,Tphylum.pvalues<0.05)
sigTphyla

# Class
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(Tclass, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0.06
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
Tmeta$Axis01 <- varespec.nmds$points[,1]
Tmeta$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(Tmeta, "Year", find_hull)
ggplot(Tmeta, aes(Axis01, Axis02), colour=Year) +
    geom_polygon(data = micro.hulls, aes(colour=Year, fill=Year), 
                 alpha = 0.1) +
    geom_point(aes(fill=Year), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ Tmeta$Year, permutations=999)
type.permanova #Significant
#Taxa driving differences
Tclass.pvalues<-apply(Tclass, 2, function(x){wilcox.test(x~Tmeta$Year)$p.value})
Tclass.pvalues.fdr<-p.adjust(Tclass.pvalues, method="fdr")
TclassPs<-data.frame(Tclass.pvalues)
sigTclasses<-subset(TclassPs,Tclass.pvalues<0.05)
sigTclasses

# Order
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(Torder, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
# NMDS Ordination
varespec.nmds <- metaMDS(varespec.bray, k=2)
varespec.nmds$stress #0
stressplot(varespec.nmds, varespec.bray)
plot(varespec.nmds, type="t", display="sites")
Tmeta$Axis01 <- varespec.nmds$points[,1]
Tmeta$Axis02 <- varespec.nmds$points[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(Tmeta, "Year", find_hull)
ggplot(Tmeta, aes(Axis01, Axis02), colour=Year) +
    geom_polygon(data = micro.hulls, aes(colour=Year, fill=Year), 
                 alpha = 0.1) +
    geom_point(aes(fill=Year), alpha=0.5, pch=21, col="black") +
    scale_size_continuous(range=c(3,10)) +
    theme_bw() +
    theme(legend.position="right")
# PERMANOVA
type.permanova<-adonis(varespec.bray ~ Tmeta$Year, permutations=999)
type.permanova #Significant
#Taxa driving differences
Torder.pvalues<-apply(Torder, 2, function(x){wilcox.test(x~Tmeta$Year)$p.value})
Torder.pvalues.fdr<-p.adjust(Torder.pvalues, method="fdr")
TorderPs<-data.frame(Torder.pvalues)
sigTorders<-subset(TorderPs,Torder.pvalues<0.05)
sigTorders

bray <- as.matrix(varespec.bray)
Tmeta <- as.matrix(Tmeta)
# write.csv(bray, "Time_Series_Bray_Curtis_Order.csv")
# write.csv(Tmeta, "Time_Series_Metadata.csv")

############################## 2007 Microbial Biomass ######################################
shapiro.test(sub$MBC)
wilcox.test(MBC ~ Class, data = sub) #NSD
aggregate(sub$MBC, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$MBC, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

################################ Enzyme Analyses ###########################################
# Not normalized. Info for supp table 2. Mean and SE
cdata <- ddply(sub, c("Class"), summarise,
               N    = length(BG),
               mean = mean(BG),
               sd   = sd(BG),
               se   = sd / sqrt(N)
)
cdata

cdata <- ddply(sub, c("Class"), summarise,
               N    = length(NAG),
               mean = mean(NAG),
               sd   = sd(NAG),
               se   = sd / sqrt(N)
)
cdata

cdata <- ddply(sub, c("Class"), summarise,
               N    = length(PHOS),
               mean = mean(PHOS),
               sd   = sd(PHOS),
               se   = sd / sqrt(N)
)
cdata

# Ratio
shapiro.test(sub$bgnag)
wilcox.test(bgnag ~ Class, data = sub)
aggregate(sub$bgnag, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$bgnag, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

shapiro.test(sub$nagphos)
wilcox.test(nagphos ~ Class, data = sub)
aggregate(sub$nagphos, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$nagphos, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

shapiro.test(sub$bgphos)
t.test(bgphos ~ Class, data = sub)
aggregate(sub$bgphos, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$bgphos, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se
# NSD for all three

# Time Series
te <- read.csv("Time_Enzymes.csv")
te$Year <- as.factor(te$Year)

shapiro.test(te$BG.NAG)
t.test(te$BG.NAG ~ te$Year)

shapiro.test(te$NAG.PHOS)
t.test(te$NAG.PHOS ~ te$Year)
# NSD for both

############################### 2007 Abiotic Factors #######################################
shapiro.test(sub$DOC)
wilcox.test(DOC ~ Class, data = sub) # p = 0.03
aggregate(sub$DOC, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$DOC, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

shapiro.test(sub$TDN)
wilcox.test(TDN ~ Class, data = sub) # p = 0.005
aggregate(sub$TDN, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$TDN, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se
ggplot(sub, aes(x=Class, y=TDN, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("TDN") +
    theme_bw()+
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))

shapiro.test(sub$Dptotal)
wilcox.test(Dptotal ~ Class, data = sub) # p = 0.01
aggregate(sub$Dptotal, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$Dptotal, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

shapiro.test(sub$PH)
t.test(PH ~ Class, data = sub) # p = 0.005
aggregate(sub$PH, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$PH, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se

shapiro.test(sub$SOIL_H2O)
wilcox.test(SOIL_H2O ~ Class, data = sub) # p = 0.88
aggregate(sub$SOIL_H2O, list(sub$Class), mean)
se <- as.data.frame(aggregate(sub$SOIL_H2O, list(sub$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sub))
se
 
shapiro.test(sub$SAND)
t.test(SAND ~ Class, data = sub) # p = 0.41
aggregate(sub$SAND, list(sub$Class), mean)

shapiro.test(sub$MeanSnow)
t.test(MeanSnow ~ Class, data = sub) # p = 0.0008
aggregate(sub$MeanSnow, list(sub$Class), mean)



################################### Soddie Data ############################################
sod <- read.csv("Soddie.csv")
tax <- read.csv("SoddieMajor.csv")
sod$bgnag <- sod$BG/sod$NAG
sod$nagphos <- sod$NAG/sod$PHOS
sod$bgphos <- sod$BG/sod$PHOS

# Specific Taxa test
wilcox.test(tax$AcidobacteriaAcidobacteriia ~ tax$site)
wilcox.test(tax$ActinobacteriaThermoleophilia ~ tax$site)
wilcox.test(tax$ActinobacteriaThermoleophiliaSolirubrobacterales ~ tax$site)
wilcox.test(tax$ChloroflexiKtedonobacteria ~ tax$site)
wilcox.test(tax$ProteobacteriaAlpha ~ tax$site)
wilcox.test(tax$AlphaRhizobiales ~ tax$site)

# Microbial Biomass
shapiro.test(sod$MBC)
t.test(MBC ~ Class, data = sod)
aggregate(sod$MBC, list(sod$Class), mean)
se <- as.data.frame(aggregate(sod$MBC, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se
ggplot(sod, aes(x=Class, y=MBC, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Microbial Biomass Carbon") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))

shapiro.test(sod$MBN)
t.test(MBN ~ Class, data = sod)
ggplot(sod, aes(x=Class, y=MBN, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Microbial Biomass Nitrogen") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))

# Soil Properties
shapiro.test(sod$ExtractC)
t.test(ExtractC ~ Class, data = sod)
se <- as.data.frame(aggregate(sod$ExtractC, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se
ggplot(sod, aes(x=Class, y=ExtractC, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Soil Extractable C") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly more soil extractable C in vascular (p=0.0008)

shapiro.test(sod$ExtractN)
t.test(ExtractN ~ Class, data = sod)
se <- as.data.frame(aggregate(sod$ExtractN, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se
ggplot(sod, aes(x=Class, y=ExtractN, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Soil Extractable N") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly more soil extractable N in vascular (p = 0.007)

shapiro.test(sod$pH)
wilcox.test(pH ~ Class, data = sod) #NSD
aggregate(sod$pH, list(sod$Class), mean)
se <- as.data.frame(aggregate(sod$pH, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se

shapiro.test(sod$Soil_H2O)
t.test(Soil_H2O ~ Class, data = sod)
se <- as.data.frame(aggregate(sod$Soil_H2O, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se
ggplot(sod, aes(x=Class, y=Soil_H2O, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Soil Moisture") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly more soil moisture in vascular (p<0.0001)

# Enzymes
# Summary of raw data for Supp Table 2
cdata <- ddply(sod, c("Class"), summarise,
               N    = length(BG),
               mean = mean(BG),
               sd   = sd(BG),
               se   = sd / sqrt(N)
)
cdata

cdata <- ddply(sod, c("Class"), summarise,
               N    = length(NAG),
               mean = mean(NAG),
               sd   = sd(NAG),
               se   = sd / sqrt(N)
)
cdata

cdata <- ddply(sod, c("Class"), summarise,
               N    = length(PHOS),
               mean = mean(PHOS),
               sd   = sd(PHOS),
               se   = sd / sqrt(N)
)
cdata

# Ratios
shapiro.test(sod$bgnag)
t.test(bgnag ~ Class, data = sod)
aggregate(sod$bgnag, list(sod$Class), mean)
se <- as.data.frame(aggregate(sod$bgnag, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se

shapiro.test(sod$nagphos)
t.test(nagphos ~ Class, data = sod) # Significant!
aggregate(sod$nagphos, list(sod$Class), mean)
se <- as.data.frame(aggregate(sod$nagphos, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se

shapiro.test(sod$bgphos)
wilcox.test(bgphos ~ Class, data = sod) # Significant!
aggregate(sod$bgphos, list(sod$Class), mean)
se <- as.data.frame(aggregate(sod$bgphos, list(sod$Class), FUN = sd))
se$se <- se$x/sqrt(NROW(sod))
se

plot(log(sod$PHOS), log(sod$BG))
plot(log(sod$PHOS), log(sod$NAG))
plot(log(sod$NAG), log(sod$BG))

## CO2 and Nfixation
shapiro.test(sod$CO2l)
wilcox.test(CO2l ~ Class, data = sod)
ggplot(sod, aes(x=Class, y=CO2l, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("CO2") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly higher CO2l in vascular (p=0.01)

shapiro.test(sod$qCO2l)
wilcox.test(qCO2l ~ Class, data = sod)

shapiro.test(sod$CO2d)
wilcox.test(CO2d ~ Class, data = sod)
ggplot(sod, aes(x=Class, y=CO2d, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("CO2") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly higher CO2d in vascular (p=0.001)

shapiro.test(sod$qCO2d)
wilcox.test(qCO2d ~ Class, data = sod)

shapiro.test(sod$ethL)
wilcox.test(ethL ~ Class, data = sod)
ggplot(sod, aes(x=Class, y=ethL, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Ethylene") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly more ethylene in vascular (p=0.04988)

shapiro.test(sod$nfixL)
wilcox.test(nfixL ~ Class, data = sod)

shapiro.test(sod$ethD)
t.test(ethD ~ Class, data = sod)
ggplot(sod, aes(x=Class, y=ethD, fill=Class)) +
    geom_jitter(alpha=0.5, size=2, pch=21) +
    geom_boxplot(alpha=0.4, outlier.shape = NA) +
    xlab(NULL) +
    ylab("Ethylene") +
    theme_bw() +
    theme(legend.position="right", axis.text.x=element_text(angle=45, hjust=1))
# Significantly more ethylene D in vascular (p=0.0045)

shapiro.test(sod$nfixD)
t.test(nfixD ~ Class, data = sod)

################################# Graphs ###################################################
# Figure 3
g <- read.csv("ChangeGraph2.csv")
ggplot(g, aes(x=Taxon, y=Difference, fill=Site)) +
    geom_bar(aes(fill=Site), position = "dodge", stat='identity') +
    coord_flip() +
    ylab("Difference in Mean Relative Abundance \n (Vascular - Moss)") +
    scale_fill_manual(values = c("gray75", "gray25"),
                      guide = guide_legend(reverse=TRUE)) +
    scale_x_discrete(limits=c("Spartobacteria",
                              "Ktedonobacteria",
                              "Alphaproteobacteria",
                              "Verrucomicrobia",
                              "Proteobacteria",
                              "Chloroflexi",
                              "Actinobacteria",
                              "Acidobacteria")) +
    theme_bw() +
    theme(legend.position = c(0.9, 0.9),
          legend.background = element_blank(),
          legend.title = element_blank(),
          axis.title.x = element_text(face="bold", size = 16, vjust = 0), 
          axis.text.x = element_text(size = 14), 
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(face="bold", size=16))


################### Full Arikaree Time Series Ordination ###################################
meta <- read.csv("Time_Metadata.csv", header = TRUE, row.names=1)
meta$Year <- as.factor(meta$Year)
order <- t(read.table("moss_tax_L4.txt", sep="\t", header=T, row.names=1))
order <- order[1:40,]
rownames(meta)==rownames(order)
# Hellinger transformation and Bray-Curtis distance
varespec.hel <- decostand(order, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")

# PCoA Figure 2
varespec.pcoa<-cmdscale(varespec.bray, k=2, eig=T)
eig2 <- eigenvals(varespec.pcoa)
eig2 / sum(eig2)
ordiplot(scores(varespec.pcoa)[,1:2], type="t", xlab="PCo1", ylab="PCo2")
meta$Axis01 <- scores(varespec.pcoa)[,1]
meta$Axis02 <- scores(varespec.pcoa)[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(meta, "Type", find_hull)
g1 <- ggplot(meta, aes(Axis01, Axis02)) +
    geom_polygon(data = micro.hulls, aes(colour=Type), fill = NA) +
    geom_point(aes(colour=factor(Type), shape = factor(Year), size = 4)) +
    scale_colour_manual(values = c("blue", "red")) +
    xlab("PC1 - Percent Variation Explained 30.45%") +
    ylab("PC2 - Percent Variation Explained 16.92%") +
    ggtitle("PCoA - PC1 vs PC2") +
    theme_bw() +
    theme(legend.position="none")

# Soddie Graph
sod <- read.csv("Soddie.csv", header = TRUE, row.names=1)
s <- t(read.table("Soddie_L4.txt", sep="\t", header=T, row.names=1))
rownames(sod)==rownames(s)
varespec.hel <- decostand(s, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
varespec.pcoa<-cmdscale(varespec.bray, k=2, eig=T)
eig2 <- eigenvals(varespec.pcoa)
eig2 / sum(eig2)
sod$Axis01 <- scores(varespec.pcoa)[,1]
sod$Axis02 <- scores(varespec.pcoa)[,2]
# Graph Ordination with Hulls
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
micro.hulls <- ddply(sod, "Class", find_hull)
g2<-ggplot(sod, aes(Axis01, Axis02)) +
    geom_polygon(data = micro.hulls, aes(colour=Class), fill = NA) +
    geom_point(aes(colour=factor(Class), size = 4)) +
    scale_colour_manual(values = c("blue", "red")) +
    xlab("PC1 - Percent Variation Explained 38.07%") +
    ylab("PC2 - Percent Variation Explained 16.44%") +
    ggtitle("PCoA - PC1 vs PC2") +
    theme_bw() +
    theme(legend.position="none")

multiplot(g2,g1,cols=2) # Figure 2
# Moss is blue on right, Vasc is red on left. Label points that changed over time.

# Control Analysis - was change in microbes due to vegetation change or just time?
d <- subset(meta, Change == "No")
vasc <- subset(d, Type == "Vascular")
vascOrder <- order[c(1,3,6,7,9,10,11,12,15,16,25,27,30:37),]
varespec.hel <- decostand(vascOrder, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
type.permanova<-adonis(varespec.bray ~ vasc$Year, permutations=999)
type.permanova #Significant
moss <- subset(d, Type == "Moss")
mossOrder <- order[c(2,4,5,18,19,20,26,28,29,38,39,40),]
varespec.hel <- decostand(mossOrder, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")
type.permanova<-adonis(varespec.bray ~ moss$Year, permutations=999)
type.permanova #Not significant - Moss plots that didn't switch to vascular did not change!

######################### Natural Logarithm Enzyme Graph ###################################
ln <- read.csv("NatLog.csv")
ln$lnBG <- log(ln$BG)
ln$lnPHOS <- log(ln$PHOS)
ln$lnNAG <- log(ln$NAG)

# SMA Regressions
m <- lmodel2(lnBG ~ lnPHOS, data=ln, nperm=1000, range.x="relative", range.y="relative")
m
m$rsquare
m$P.param

m1 <- lmodel2(lnNAG ~ lnPHOS, data=ln, nperm=1000, range.x="relative", range.y="relative")
m1
m1$rsquare
m1$P.param

m2 <- lmodel2(lnBG ~ lnNAG, data=ln, nperm=1000, range.x="relative", range.y="relative")
m2
m2$rsquare
m2$P.param

g1 <- ggplot(ln, aes(x=lnPHOS, y=lnBG)) +
    geom_point(aes(colour=factor(Class), shape = factor(Site)), size = 4) +
    scale_colour_manual(values = c("blue", "red"), 
                        guide = guide_legend(reverse=TRUE,
                                             override.aes=list(shape=15,
                                                               shape=15))) +
    geom_abline(slope = 0.9384502, intercept = -0.20459725) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    xlab("ln(AP)") +
    ylab("ln(BG)") +
    ylim(0,8) + 
    xlim(0,8) +
    theme_bw() +
    theme(legend.position=c(0.25,0.83), legend.title=element_blank(), 
          legend.key = element_blank(), legend.background = element_blank())

g2 <- ggplot(ln, aes(x=lnPHOS, y=lnNAG)) +
    geom_point(aes(colour=factor(Class), shape = factor(Site), size = 4)) +
    scale_colour_manual(values = c("blue", "red")) +
    geom_abline(slope = 1.197110, intercept = -3.104195) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    xlab("ln(AP)") +
    ylab("ln(NAG)") +
    ylim(0,8) + 
    xlim(0,8) +
    theme_bw() +
    theme(legend.position="null")

g3 <- ggplot(ln, aes(x=lnNAG, y=lnBG)) +
    geom_point(aes(colour=factor(Class), shape = factor(Site)), size = 4) +
    scale_colour_manual(values = c("blue", "red")) +
    geom_abline(slope = 0.7839300, intercept = 2.228875) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    xlab("ln(NAG)") +
    ylab("ln(BG)") +
    ylim(0,8) + 
    xlim(0,8) +
    theme_bw() +
    theme(legend.position="null")

multiplot(g1,g2,g3,cols=3) # Figure 4

###################### Test for any differences in sequencing depth ########################
d <- read.csv("SequenceDepth.csv")
d$Year <- as.factor(d$Year)
shapiro.test(d$Sequences)
hist(d$Sequences)
t.test(d$Sequences ~ d$Year) # NSD
t.test(d$Sequences ~ d$Type) # NSD
