# Analysis of root microscopy data across the Green Lakes Valley, Colorado
# Samples collected from grid of plots established by Andrew King in 2007
# Now published in Fungal Ecology
# 1. Test AMF and DSE colonization for spatial autocorrelation
# 2. Correlations between AMF and DSE
# 3. Analyze AMF and DSE colonization by plant species and functional groups
# 4. Find best combination of predictors of AMF and DSE colonization
# 5. Plant phylogenetic analysis
# 6. Analyze fungal community composition from sequence data

######################################## Setup ############################################
library(ggplot2)
library(lmodel2)
library(Matrix)
library(lme4)
library(leaps)
library(bestglm)
library(AICcmodavg)
library(lattice)
library(survival)
library(Formula)
library(Hmisc)
library(PMCMR)
library(reshape)
library(corrplot)
library(MASS)
library(grid)
library(vcd)
library(car)
library(nlme)
library(pwr)
library(ape)
library(picante)
library(phytools)
library(rsq)
library(optimx)
library(modEvA)
library(permute)
library(vegan)
library(plyr)
library(maps)
library(indicspecies)
library(cluster)
library(rcompanion)
source("~/Desktop/Functions/Summary.R")
source("~/Desktop/Functions/logisticPseudoR2s.R")
setwd("~/Desktop/CU/2Research/RootSamples")
d <- read.csv("AK_Root_Colonization.csv")
d$Plot <- as.factor(d$Plot)

# N:P Ratio
plot(d$N.P, d$DSE)
m <- lm(d$DSE ~ d$N.P)
summary(m)
plot(d$N.P, d$AM)
m <- lm(d$AM ~ d$N.P)
summary(m)
plot(d$TIN, d$N.P)
plot(d$TDIP, d$N.P)

# Data transformations
# Ranelli et al. 2015 logit transformed fungal data and nat log transformed some predictors
# I will logit transform AMF and DSE colonization and log transform plant density, TIN, and DIP
d$AMlogit <- logit(d$AM)
d$DSElogit <- logit(d$DSE)
d$LogDensity <- log(d$Density)
d$LogTDIP <- log(d$TDIP)
d$LogTIN <- log(d$TIN)
d$LogN.P <- log(d$N.P)

# Check correlations between Elevation, Snowpack, Density, TDP, TIN (looks okay)
env <- d[,c(27,87,23,34,91)]
env <- subset(env, TDP != "NA")
M <- cor(env)
corrplot(M, method = "number", type = "lower")

############################# Spatial autocorrelation #####################################
root.dists <- as.matrix(dist(cbind(d$Easting, d$Northing)))
root.dists.inv <- 1/root.dists
root.dists.inv[1:5, 1:5]
root.dists.inv <- replace(root.dists.inv, root.dists.inv == "Inf", 0)
Moran.I(d$AM, root.dists.inv)
Moran.I(d$DSE, root.dists.inv)

############################# AMF and DSE Correlations ####################################
# Full dataset
cor.test(d$DSE, d$AM, alternative = "two.sided", method = "spearman") # NSD
plot(d$DSE, d$AM) # Looks like negative exponential

# Estimate lambda
attach(d)
data <- data.frame(AM, DSE) 
nll <- function (theta, data) {
  lambda <- theta[1]
  y <- data[,1]
  d <- data[,2]
  mu <- exp(-lambda*d)
  sigma <- exp(theta[2])
  -sum(dnorm(y, mu, sigma, log = TRUE))
}
theta_init <- c(lambda = 25, log_sigma = 0) 
res <- optim(theta_init, nll, data = data)
res
MLE_optim <- c(res$par[1], res$par[2]) 
MLE_optim
Fitted <- exp(-MLE_optim[1]*data[,2]) 
Residuals <- data[,1] - Fitted
plot(Fitted, Residuals)
sorted <- data[order(data[,2]),]
fitted <- exp(-MLE_optim[1]*sorted[,2])
plot(sorted[,2], sorted[,1], xlab = "DSE", ylab = "AMF")
lines(sorted[,2], fitted, col = "Blue", lwd = 2.5)

# Use nls
fit_nls = nls(AM ~ exp(-b*DSE), 
              start=list(b=19),
              lower=list(b=19),
              upper=list(b=19),
              algorithm = "port")
summary(fit_nls) # NS
detach(d)

# Only look at roots that had both
both <- subset(d, AM > 0 & DSE > 0)
cor.test(both$DSE, both$AM, alternative = "two.sided", method = "spearman") # NSD

attach(both)
data <- data.frame(AM, DSE) 
nll <- function (theta, data) {
  lambda <- theta[1]
  y <- data[,1]
  d <- data[,2]
  mu <- exp(-lambda*d)
  sigma <- exp(theta[2])
  -sum(dnorm(y, mu, sigma, log = TRUE))
}
theta_init <- c(lambda = 25, log_sigma = 0) 
res <- optim(theta_init, nll, data = data)
res
MLE_optim <- c(res$par[1], res$par[2]) 
MLE_optim
Fitted <- exp(-MLE_optim[1]*data[,2]) 
Residuals <- data[,1] - Fitted
plot(Fitted, Residuals)
sorted <- data[order(data[,2]),]
fitted <- exp(-MLE_optim[1]*sorted[,2])
plot(sorted[,2], sorted[,1], xlab = "DSE", ylab = "AMF")
lines(sorted[,2], fitted, col = "Blue", lwd = 2.5)

fit_nls = nls(AM ~ exp(-b*DSE), 
              start=list(b=19.902954),
              lower=list(b=19.902954),
              upper=list(b=19.902954),
              algorithm = "port")
summary(fit_nls) # NS
detach(both)

# Figure 6
ggplot(both,aes(DSE, AM)) +
  geom_point(pch = 1, size = 3, alpha = 0.75) +
  xlab("% DSE Colonization") +
  ylab("% AMF Colonizaiton") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold",size = 14), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(face = "bold",size=14),
        axis.text.y = element_text(size = 12))

######################### Plant Species and Functional Group ##############################
# Analysis with species with 3 or more samples
toBeRemoved<-which(d$Species=="AGRVAR"|d$Species=="AQUCOE"|d$Species=="CARHET"|d$Species=="CARPER"|d$Species=="CARRUP"|d$Species=="CARSCO"|d$Species=="CARSPP"|d$Species=="CERARV"|d$Species=="LLOSER"|d$Species=="LUZSPI"|d$Species=="MERLAN"|d$Species=="OREALP"|d$Species=="SIBPRO"|d$Species=="UNKOPP")
n3 <- droplevels(d[-toBeRemoved, ])
# write.csv(n3, "n3.csv")
# write.csv(n3, "n3orig.csv")
sum <- summarySE(n3, measurevar = "AM", groupvars = c("Species"))
sum <- summarySE(n3, measurevar = "DSE", groupvars = c("Species"))

# Test for Homoscedasticity
leveneTest(n3$AM ~ n3$Species) # p < 0.05, not homogeneous
leveneTest(n3$DSE ~ n3$Species) # p < 0.05, not homogeneous

# Species
kruskal.test(n3$AM ~ n3$Species) # Significant
kruskal.test(n3$DSE ~ n3$Species) # Significant

# Pairwise comparisons
posthoc.kruskal.nemenyi.test(n3$AM ~ n3$Species)
posthoc.kruskal.nemenyi.test(n3$DSE ~ n3$Species)

# Combined Flipped, for Tree
g <- read.csv("n3.csv")
toBeRemoved<-which(g$Species=="POASPP")
g <- droplevels(g[-toBeRemoved, ])
fungi_names <- c(`AMF` = "a) AMF", `DSE` = "b) DSE")
ggplot(g,aes(reorder(Species2, Tree, median), Colonization, colour = Group)) + 
  scale_x_discrete(limits=c("EriSim","SenFre","AntMed","CirSco","AngGra","BesAlp","MinObt","SilAca","SteUmb","OxyDig","TriDas","GeuRos","CarAlb","CarPha","CarPyr","KobMyo","DesCes","ElyScr","FesBra","TriSpi")) +
  geom_boxplot() +
  facet_wrap(~Fungus, labeller = as_labeller(fungi_names)) +
  xlab(NULL) +
  ylab("Percent Colonization") +
  theme_bw() +
  theme(legend.position = c(0.42,0.8),
        legend.background = element_rect(colour = "black"),
        axis.title.x = element_text(face="bold",size = 16, vjust = 0), 
        axis.text.x = element_text(size = 12,hjust = 0.5), 
        axis.text.y = element_text(size = 14,vjust=0.5,hjust=0), 
        axis.title.y = element_text(face="bold", size=16),
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size = 16)) +
  coord_flip() +
  ylim(0,75)

# Forb/Grass/Sedge (Functional Group)
# Only 4 samples for N-fixer so remove
# Only 2 samples for rush so remove
toBeRemoved<-which(d$Species=="LUZSPI"|d$Species=="TRIDAS")
fg <- droplevels(d[-toBeRemoved, ])
write.csv(fg, "fg.csv")

# Test for Normal distribution
shapiro.test(fg$AM) # p < 0.05, not normal
shapiro.test(fg$DSE) # p < 0.05, not normal

leveneTest(fg$AM ~ fg$Type) # p < 0.05, not homogeneous 
leveneTest(fg$DSE ~ fg$Type) # p < 0.05, not homogeneous

kruskal.test(fg$AM ~ fg$Type) # Significant
kruskal.test(fg$DSE ~ fg$Type) # Significant

# Pairwise comparisons
posthoc.kruskal.nemenyi.test(fg$AM ~ fg$Type)
# Label graph a, ab, b
posthoc.kruskal.nemenyi.test(fg$DSE ~ fg$Type)
# Label graph a, a, b

# Combined
g2 <- read.csv("fg.csv")
ggplot(g2,aes(Type, Colonization)) + 
  geom_boxplot() +
  facet_wrap(~Fungus, labeller = as_labeller(fungi_names)) +
  xlab("Functional Group") +
  ylab("Percent Colonization") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size = 16, vjust = 0), 
        axis.text.x = element_text(angle = 90,size = 12), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold", size=16),
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size = 16)) +
  ylim(0,75)

# Functional Group where species average is data point
# Forb vs. Grass vs. Sedge
toBeRemoved<-which(d$Species=="LUZSPI"|d$Species=="TRIDAS")
fgs <- droplevels(d[-toBeRemoved, ])

AMsum <- summarySE(fgs,measurevar = "AM", groupvars = c("Species", "Type"))
table(AMsum$Type)
leveneTest(AMsum$AM ~ AMsum$Type) # Homogeneous
m <- aov(AMsum$AM ~ AMsum$Type)
summary(m) # NSD
shapiro.test(m$residuals) # Not normal
kruskal.test(AMsum$AM ~ AMsum$Type) # Marginal
posthoc.kruskal.nemenyi.test(AMsum$AM ~ AMsum$Type) # NSD

DSsum <- summarySE(fgs,measurevar = "DSE", groupvars = c("Species", "Type"))
table(DSsum$Type)
leveneTest(DSsum$DSE ~ DSsum$Type) # Homogeneous
m <- aov(DSsum$DSE ~ DSsum$Type)
summary(m)
shapiro.test(m$residuals) # Not normal
kruskal.test(DSsum$DSE ~ DSsum$Type) # Significant
posthoc.kruskal.nemenyi.test(DSsum$DSE ~ DSsum$Type) # AB, A, B, sedge lowest

# Figure 4
func <- read.csv("Species_Means.csv")
ggplot(func, aes(Type, Colonization)) +
  geom_boxplot() +
  facet_wrap(~Fungus, labeller = as_labeller(fungi_names)) +
  xlab("Functional Group") +
  ylab("Percent Colonization") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size = 16, vjust = 0), 
        axis.text.x = element_text(angle = 90,size = 12), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold", size=16),
        strip.background = element_blank(),
        strip.text = element_text(face="bold",size = 16))



###################### Multivariate Best Combination of Predictors ########################
# These models can't handle, NA's, so make a new dataframe without NA's
d <- subset(d, TIN != "NA")
d <- subset(d, TDIP != "NA")
d <- subset(d, LogTDIP != "-Inf")
d <- subset(d, LogN.P != "-Inf")

# Logit Y, Species, Elevation, Mean.Snow, LogDensity, LogTDIP, LogTIN, LogN.P
env <- d[,c(27,87,102,103,104,105)]
M <- cor(env)
corrplot(M, method = "number", type = "lower")
X <- as.data.frame(d[,c(15,27,87,102,103,104,105)])
y <- d$AMlogit
Xy <- as.data.frame(cbind(X,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel
bestmodel$BestModels
nullmodel <- lm(AMlogit ~ 1, data = d)
bestmodel <- lm(AMlogit ~ Species + TIN, data = d)
summary(bestmodel)
# plot(bestmodel) # 27, 83, 139, 142
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
# Partial R2 for top 5
bestmodel <- lm(AMlogit ~ Species + LogTIN + Elevation + Snowpack, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(AMlogit ~ Species + LogTIN + Elevation + Snowpack + LogTDIP, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(AMlogit ~ Species + LogTIN + Elevation, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(AMlogit ~ Species + LogTIN + Elevation + LogTDIP, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(AMlogit ~ Species + LogTIN, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)

X <- as.data.frame(d[,c(15,27,87,102,103,104,105)])
y <- d$DSElogit
Xy <- as.data.frame(cbind(X,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel
bestmodel$BestModels
nullmodel <- lm(DSElogit ~ 1, data = d)
bestmodel <- lm(DSElogit ~ Mean.Snow + Elevation, data = d)
summary(bestmodel)
# plot(bestmodel) # 15, 18, 138
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
bestmodel <- lm(DSElogit ~ Snowpack + LogTIN, data = d)
summary(bestmodel)
# Partial R2 for top 5
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(DSElogit ~ Snowpack + LogTIN + LogDensity, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(DSElogit ~ Snowpack + LogDensity, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(DSElogit ~ Snowpack, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)
bestmodel <- lm(DSElogit ~ Snowpack + Elevation, data = d)
summary(bestmodel)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)

# Figure 5 insets
ggplot(d,aes(TDIP, AM)) +
  geom_point(pch = 1, size = 3, alpha = 0.75) +
  xlab("TDIP (mg/L)") +
  ylab("% AMF") +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold",size = 14), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(face = "bold",size=14),
        axis.text.y = element_text(size = 12))

ggplot(d,aes(Mean.Snow, DSE)) +
  geom_point(pch = 1, size = 3, alpha = 0.75) +
  xlab("Snowpack (cm)") +
  ylab("% DSE") +
  scale_x_continuous(breaks = c(100, 200, 300)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold",size = 14), 
        axis.text.x = element_text(size = 12,hjust = 0.65),
        axis.title.y = element_text(face = "bold",size=14),
        axis.text.y = element_text(size = 12))

############################ Phylogenetic Analyses #########################################
# with some code and help from Marko Spasojevic
tre=read.newick(text = "((((((((((((((((((((((((((((((((((((((Erigeronsimplex:22.8719)Erigeron:22.8719,(Seneciofremontii:22.8719)Senecio:22.8719,(Antennariamedia:22.8719)Antennaria:22.8719,(Cirsiumscopulorum:22.8719)Cirsium:22.8719)Asteraceae:3.346060):7.104420):16.236000):6.015739):0.292069):4.165166):0.119934):0.640582)Asterales:11.753099,((((((((Angelicagrayi:35.7241)Angelica:35.7241,(Oreoxisalpina:35.7241)Oreoxis:35.7241)Apiaceae:1.941560):2.849080):0.340050):4.487350):1.314090):9.328050)Apiales:3.708482):0.900124):0.056880):0.252052):0.788782):3.574670)Campanulidae:1.882210,(((((((Mertensialanceolata:27.8976)Mertensia:27.8976)Boraginaceae:19.203790,((((((((((((((Besseyaalpina:18.7289)Besseya:18.7289)Scrophulariaceae:3.509661):0.094686):0.332342):0.621584):0.129732):0.182246):0.432809):0.696615):12.782970):5.153970)Lamiales:11.432543):0.542387):1.629560):0.720309):3.417357):17.813603):4.580345)Lamiidae:1.340875):4.737980)Asteridae:8.908438):0.179371,(((((((((Cerastiumarvense:25.3025)Cerastium:25.3025,(Minuartiaobtusiloba:25.3025)Minuartia:25.3025,(Sileneacaulis:25.3025)Silene:25.3025,(Stellariaumbellata:25.3025)Stellaria:25.3025)Caryophyllaceae:6.366525):4.028800):16.488400):5.809780):3.749530):2.412990,(((((Oxyriadigyna:27.0072)Oxyria:27.0072)Polygonaceae:5.087963):24.304300):2.270850):3.783420)Caryophyllales:25.700064):1.536275)Pentapetalae:1.460940)Superasteridae:0.986915):0.020375,((((((((((Trifoliumdasyphyllum:34.5713)Trifolium:34.5713)Fabaceae:0.276658):2.378927):0.250678)Fabales:34.946600,((((Geumrossii:41.147)Geum:41.147,(Sibbaldiaprocumbens:41.147)Sibbaldia:41.147)Rosaceae:6.983035)Rosales:16.658900):1.059410):6.263220)Fabidae:0.420986):3.855140)Rosidae:0.138911)Superrosidae:1.491856)Gunneridae:17.734960,((((((((Aquilegiacoerulea:52.7228)Aquilegia:52.7228)Ranunculaceae:18.954800):3.572080):2.332450):2.577480):0.137456)Ranunculales:3.072206):0.808344)Eudicotyledoneae:43.727816):7.654960,(((((((((((((((Carexscopulorum:18.1409,Carexalbonigra:18.1409,Carexheteroneura:18.1409,Carexperglobosa:18.1409,Carexphaeocephala:18.1409,Carexpyrenaica:18.1409,Carexrupestris:18.1409)Carex:18.1409,(Kobresiamyosuroides:18.1409)Kobresia:18.1409)Cyperaceae:12.397796,((Luzulaspicata:16.2789)Luzula:16.2789)Juncaceae:16.121726):8.029155):9.161151):11.650487):6.700159,(((((((Agrostisvariabilis:23.8819)Agrostis:23.8819,(Deschampsiacespitosa:23.8819)Deschampsia:23.8819,(Elymusscriberneri:23.8819)Elymus:23.8819,(Festucabrachyphylla:23.8819)Festuca:23.8819,(Trisetumspicatum:23.8819)Trisetum:23.8819)Poaceae:1.789542):4.309137):7.809344):6.449827):13.917930):2.180905):1.372227):4.576783)Poales:33.002000)Commelinidae:11.693860,((((((((Lloydiaserotina:29.2241)Lloydia:29.2241)Liliaceae:20.642664):1.497000):25.871794):6.208706):16.804775)Liliales:5.239854):0.152371):2.668590):19.828273)Petrosaviidae:0.979827)Nartheciidae:12.924900)Monocotyledoneae:17.016266):2.257320):3.495264)Mesangiospermae:23.596800):16.908091):8.728959)Angiospermae:108.965000)Spermatophyta:38.467865):10.085114):4.856365):5.208675):11.765386):17.675263;")

tre
tree<-collapse.singles(tre)
tree
plot(tree)
plot(tree,type="fan",cex = 0.75)
plot(tree, direction = "upwards", cex = 0.9)
p.dist <- cophenetic.phylo(tree)
# Prune tree for sp sampled 3 times
phy<-collapse.singles(tre)
length(phy$tip.label)
phy$Nnode
comm<-matrix(nrow=1,ncol=20)
colnames(comm) <- c("Angelicagrayi","Antennariamedia","Besseyaalpina","Carexalbonigra","Carexphaeocephala","Carexpyrenaica","Cirsiumscopulorum","Deschampsiacespitosa","Elymusscriberneri","Erigeronsimplex","Festucabrachyphylla","Geumrossii","Kobresiamyosuroides","Minuartiaobtusiloba","Oxyriadigyna","Seneciofremontii","Sileneacaulis","Stellariaumbellata","Trifoliumdasyphyllum","Trisetumspicatum")
prunedphy <- prune.sample(comm, phy)
plot(prunedphy)

# make tree dichotamous
prunedphy2 <- multi2di(prunedphy)
plot(prunedphy2)

# Say mean DSE or AMF is a trait
DSEtrait <- c(10.3,0.5,0.25,0,1.7,1.6,0.5,27.5,0.4,2.3,3.3,8.7,0.6,0.3,0.1,1,1.4,4.3,8,11)
AMFtrait <- c(11,5.5,24,20.6,18,42.3,0.2,6.8,0.1,4.2,0,8.3,14,23,2.8,1.1,9.3,1.8,0.9,2.5)
traits <- data.frame(DSEtrait,AMFtrait)
rownames(traits) <- c("Erigeronsimplex","Seneciofremontii","Antennariamedia","Cirsiumscopulorum","Angelicagrayi","Besseyaalpina","Minuartiaobtusiloba","Sileneacaulis","Stellariaumbellata","Oxyriadigyna","Trifoliumdasyphyllum","Geumrossii","Carexalbonigra","Carexphaeocephala","Carexpyrenaica","Kobresiamyosuroides","Deschampsiacespitosa","Elymusscriberneri","Festucabrachyphylla","Trisetumspicatum")
multiPhylosignal(traits,prunedphy2)

############################ Fungal Composition Analysis ###################################
# Genera Richness (from ITS sequences) info for Table 1
aggregate(AM_Rich~Species, FUN=mean, data = d)
richmin <- aggregate(AM_Rich~Species, FUN=min, data = d)
richmax <- aggregate(AM_Rich~Species, FUN=max, data = d)
range <- data.frame(richmin, richmax)
range
aggregate(DSE_Rich~Species, FUN=mean, data = d)
richmin <- aggregate(DSE_Rich~Species, FUN=min, data = d)
richmax <- aggregate(DSE_Rich~Species, FUN=max, data = d)
range <- data.frame(richmin, richmax)
range

# Analysis of AMF and DSE community composition, genus level
# ITS sequence samples from 2016 roots, rarefied at 5238
# n = 136 that passed rarefaction cutoff and had environmental data
# Predictors: Elevation, Snowpack (mean 1997-2015), TDIP (2007), TIN (2015), plant density (2015)
# Also test plant species and functional group

### AMF
setwd("~/Desktop/CU/2Research/RootSamples")
meta <- read.csv("AMF_Comp_Meta_NoNA.csv")
AMF_genus <- read.csv("AMF_Genera_Comp_NoNA_names.csv")
AMF_genus$Code == meta$Code
AMF_genus <- AMF_genus[,-1]

# Rank Abundance Curve
spe.tot <- apply(AMF_genus, 2, sum)
sort(spe.tot)
plot(radfit(spe.tot)) # Lognormal

# Hellinger, Bray Curtis Matrix
varespec.hel <- decostand(AMF_genus, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")

# Principle Coordinates Analysis
varespec.pcoa<-cmdscale(varespec.bray, k=nrow(AMF_genus)-1, eig=T)
eig2 <- eigenvals(varespec.pcoa)
eig2 / sum(eig2)
meta$fungAxis01 <- scores(varespec.pcoa)[,1]
meta$fungAxis02 <- scores(varespec.pcoa)[,2]
env <- scale(meta[,c(7:12)])
ef<-envfit(varespec.pcoa, env, permu=999)
ef
vec.sp.df<-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.sp.df$variables<-rownames(vec.sp.df)
vec.sp.df <- vec.sp.df[c(1,2,5,6),]
# Figure 7A
ggplot(meta, aes(fungAxis01, fungAxis02)) +
  geom_point(size = 2) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=Dim1,y=0,yend=Dim2),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               inherit.aes = FALSE) + 
  geom_text(data=vec.sp.df,aes(x=Dim1,y=Dim2,label=variables),size=4) +
  xlab("PC1 - Percent Variation Explained 22.86%") +
  ylab("PC2 - Percent Variation Explained 13.90%") +
  ggtitle("AMF PCoA - PC1 vs PC2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

m<-adonis(varespec.bray ~ meta$Species, permutations=999)
m # Not significant
m<-betadisper(varespec.bray, meta$Species)
anova(m) # Significant

# Group - n = 1 luzula as grass
m<-adonis(varespec.bray ~ meta$Group, permutations=999)
m # Not significant
m<-betadisper(varespec.bray, meta$Group)
anova(m) # Not Significant

###DSE
setwd("~/Desktop/CU/2Research/RootSamples")
meta <- read.csv("DSE_Comp_Meta_NoNA.csv")
DSE_genus <- read.csv("DSE_Genera_Comp_NoNA_names.csv")
DSE_genus$Code == meta$Code
DSE_genus <- DSE_genus[,-1]

# Rank Abundance Curve
spe.tot <- apply(DSE_genus, 2, sum)
sort(spe.tot)
plot(radfit(spe.tot)) # Lognormal

# Hellinger, Bray Curtis Matrix
varespec.hel <- decostand(DSE_genus, method="hellinger")
varespec.bray <- vegdist(varespec.hel, method="bray")

# Principle Coordinates Analysis
varespec.pcoa<-cmdscale(varespec.bray, k=nrow(DSE_genus)-1, eig=T)
eig2 <- eigenvals(varespec.pcoa)
eig2 / sum(eig2)
meta$fungAxis01 <- scores(varespec.pcoa)[,1]
meta$fungAxis02 <- scores(varespec.pcoa)[,2]
env <- scale(meta[,c(7:12)])
ef<-envfit(varespec.pcoa, env, permu=999)
ef
vec.sp.df<-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
vec.sp.df$variables<-rownames(vec.sp.df)
vec.sp.df <- vec.sp.df[c(3,5,6),]
# Figure 7B
ggplot(meta, aes(fungAxis01, fungAxis02)) +
  geom_point(size = 2) +
  geom_segment(data=vec.sp.df,aes(x=0,xend=Dim1,y=0,yend=Dim2),
               arrow = arrow(length = unit(0.5, "cm")),colour="blue",
               inherit.aes = FALSE) + 
  geom_text(data=vec.sp.df,aes(x=Dim1,y=Dim2,label=variables),size=4) +
  xlab("PC1 - Percent Variation Explained 44.17%") +
  ylab("PC2 - Percent Variation Explained 19.03%") +
  ggtitle("DSE PCoA - PC1 vs PC2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

m<-adonis(varespec.bray ~ meta$Species, permutations=999)
m # Significant
m<-betadisper(varespec.bray, meta$Species)
anova(m) # Not Significant

# Group with n = 1 luzula as grass
m<-adonis(varespec.bray ~ meta$Group, permutations=999)
m # Significant
m<-betadisper(varespec.bray, meta$Group)
anova(m) # Not Significant
