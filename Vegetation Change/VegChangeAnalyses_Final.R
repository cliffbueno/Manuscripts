# Analysis of remote sensing data, ground truth data, and summer climate data, Niwot Ridge Colorado
# By Cliff Bueno de Mesquita Fall 2015 - Spring 2018
# Original remote sensing by Luke Tillmann 2014
# Ground truthing by Cliff Bueno de Mesquita, Connor Bernard, Katherine Rosemond
# Now published in Arctic, Antarctic, and Alpine Research

################################### Setup ##################################################
library(leaps)
library(bestglm)
library(AICcmodavg)
library(corrplot)
library(randomForest)
library(miscTools)
library(ggplot2)
library(MASS)
library(spData)
library(sp)
library(spgwr)
library(boot)
library(modEvA)
library(VSURF)
library(maptools)
library(robustbase)
library(Rcpp)
library(spdep)
library(Matrix)
library(GWmodel)
library(TTR)
library(quantmod)
library(tseries)
library(fracdiff)
library(timeDate)
library(forecast)
library(ggplot2)
library(nlme)
source("~/Desktop/Functions/logisticPseudoR2s.R")
setwd("~/Desktop/CU/2Research/VegChange")
vc <- read.csv("VegChangeAll.csv")

# What variables are too correlated with elevation? None!
plot(vc$Elevation, vc$Aspect)
plot(vc$Elevation, vc$Slope)
plot(vc$Elevation, vc$Mean)
plot(vc$Elevation, vc$CV)
plot(vc$Elevation, vc$Trend)
plot(vc$Elevation, vc$Solar)
ggplot(vc, aes(Elevation, Mean)) +
  geom_smooth(method = loess)

# Make continuous cover variable
cover <- read.csv("Elev_Cover.csv")
tcover <- cover[1:5,]
treemodel <- lm(tcover$Tree~poly(tcover$Elevation,3,raw=TRUE))
summary(treemodel)
scover <- cover[1:7,]
shrubmodel <- lm(scover$Shrub~poly(scover$Elevation,4,raw=TRUE))
summary(shrubmodel)
tundramodel <- lm(cover$Tundra~poly(cover$Elevation,4,raw=TRUE))
summary(tundramodel)

third_order <- function(newdist, model) {
    coefs <- coef(model)
    res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + 
        (coefs[4] * newdist^3)
    return(res)
}

fourth_order <- function(newdist, model) {
    coefs <- coef(model)
    res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + 
        (coefs[4] * newdist^3) + (coefs[5] * newdist^4)
    return(res)
}

vc$Tree_Cover <- third_order(vc$Elevation, treemodel)
for (i in 1:1532) {
    if (vc$Tree_Cover[i] < 0) {
        vc$Tree_Cover[i] <- 0
    }
    if (vc$Elevation[i] > 3472) {
        vc$Tree_Cover[i] <- 0
    }
}
vc$Shrub_Cover <- fourth_order(vc$Elevation, shrubmodel)
for (i in 1:1532) {
    if (vc$Shrub_Cover[i] < 0) {
        vc$Shrub_Cover[i] <- 0
    }
    if (vc$Elevation[i] > 3685) {
        vc$Shrub_Cover[i] <- 0
    }
}
vc$Tundra_Cover <- fourth_order(vc$Elevation, tundramodel)
for (i in 1:1532) {
    if (vc$Tundra_Cover[i] < 0) {
        vc$Tundra_Cover[i] <- 0
    }
}
plot(vc$Elevation, vc$Tree_Cover)
plot(vc$Elevation, vc$Shrub_Cover)
plot(vc$Elevation, vc$Tundra_Cover)

# Supplementary Figure A1
ggplot(vc, aes(x = Elevation, y = Tree_Cover)) +
  geom_line(aes(x=Elevation,y=Tree_Cover,colour="darkgreen")) +
  geom_line(aes(x=Elevation,y=Shrub_Cover,colour="lightgreen")) +
  geom_line(aes(x=Elevation,y=Tundra_Cover,colour="orange")) +
  ylim(0,100) +
  xlab("Elevation (m)") +
  ylab("Percent Cover") +
  scale_colour_manual(name = "Veg. Type", 
                      values = c("darkgreen","lightgreen","orange"),
                      labels = c("Tree","Shrub","Tundra")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"))

# covers <- as.data.frame(cbind(vc$X, vc$Y, vc$Tree_Cover, vc$Shrub_Cover, vc$Tundra_Cover))
# write.csv(covers, file = "Covers.csv")

# Subsets for tundra, shrub, open forest, absent in 1972
ta <- subset(vc, Cover_1972 != "T")
sa <- subset(vc, Cover_1972 != "S")
oa <- subset(vc, Cover_1972 != "O")
loa <- subset(oa, Elevation < 3600)
lsa <- subset(sa, Elevation < 3760)
hta <- subset(ta, Elevation > 3550)



############################ Elevational Range #############################################
# Test for change in max and min, and in 95th percentile for each
# All Forest
f72 <- subset(vc, Cover_1972 == "O" | Cover_1972 == "C")
f08 <- subset(vc, Cover_2008 == "O" | Cover_2008 == "C")
min(f72$Elevation)
max(f72$Elevation)
quantile(f72$Elevation, c(0.05, 0.95))
min(f08$Elevation)
max(f08$Elevation)
quantile(f08$Elevation, c(0.05, 0.95))

# Closed Canopy Forest
cf72 <- subset(vc, Cover_1972 == "C")
cf08 <- subset(vc, Cover_2008 == "C")
min(cf72$Elevation)
max(cf72$Elevation)
quantile(cf72$Elevation, c(0.05, 0.95))
min(cf08$Elevation)
max(cf08$Elevation)
quantile(cf08$Elevation, c(0.05, 0.95))

# Open Forest
of72 <- subset(vc, Cover_1972 == "O")
of08 <- subset(vc, Cover_2008 == "O")
min(of72$Elevation)
max(of72$Elevation)
quantile(of72$Elevation, c(0.05, 0.95))
min(of08$Elevation)
max(of08$Elevation)
quantile(of08$Elevation, c(0.05, 0.95))

# Shrub
s72 <- subset(vc, Cover_1972 == "S")
s08 <- subset(vc, Cover_2008 == "S")
min(s72$Elevation)
max(s72$Elevation)
quantile(s72$Elevation, c(0.05, 0.95))
min(s08$Elevation)
max(s08$Elevation)
quantile(s08$Elevation, c(0.05, 0.95))

# Tundra
t72 <- subset(vc, Cover_1972 == "T")
t08 <- subset(vc, Cover_2008 == "T")
min(t72$Elevation)
max(t72$Elevation)
quantile(t72$Elevation, c(0.05, 0.95))
min(t08$Elevation)
max(t08$Elevation)
quantile(t08$Elevation, c(0.05, 0.95))



################################### Best GLMs #############################################
# To test for best combo of fine-scale predictor variables of change
# First look at correlations
env <- vc[,c(3:8,27:29)]
M <- cor(env)
corrplot(M, method = "number", type = "lower")

# Tree below 3600m    
X <- as.data.frame(scale(loa[,c(3:8,27)]))
y <- loa$Oexpand
Xy <- as.data.frame(cbind(X,y))
bestLOE <- bestglm(Xy, IC = "AIC", family = binomial)
bestLOE # Cover, Elevation, Solar
bestLOE$BestModels
bestLOE <- glm(Oexpand ~ Tree_Cover + Elevation + Solar, family = binomial, data = loa)
summary(bestLOE)
logisticPseudoR2s(bestLOE)
Dsquared(bestLOE, adjust = TRUE)
bestLOECV<-cv.glm(data=loa,glmfit=bestLOE,K=10)
bestLOECV$delta
bestLOEAUC<-AUC(model=bestLOE)
bestLOEAUC$AUC
# Null AIC 318.88, AIC 278.27
min(loa$Tree_Cover)
max(loa$Tree_Cover)
min(loa$Elevation)
max(loa$Elevation)
min(loa$Solar)
max(loa$Solar)

# Shrub below 3760m
lsa <- subset(sa, Elevation < 3760)
X <- as.data.frame(scale(lsa[,c(3:8,28)]))
y <- lsa$Sexpand
Xy <- as.data.frame(cbind(X,y))
bestLSE <- bestglm(Xy, IC = "AIC", family = binomial)
bestLSE # Cover, Elevation, Solar, Trend
bestLSE$BestModels
bestLSE <- glm(Sexpand ~ Shrub_Cover + Elevation + Solar + Trend, family = binomial, data = lsa)
summary(bestLSE)
logisticPseudoR2s(bestLSE)
Dsquared(bestLSE, adjust = TRUE)
bestLSECV<-cv.glm(data=lsa,glmfit=bestLSE,K=10)
bestLSECV$delta
bestLSEAUC<-AUC(model=bestLSE)
bestLSEAUC$AUC
# Null AIC 352.63, AIC 321.65
min(lsa$Shrub_Cover)
max(lsa$Shrub_Cover)
min(lsa$Elevation)
max(lsa$Elevation)
min(lsa$Solar)
max(lsa$Solar)
min(lsa$Trend)
max(lsa$Trend)

# Tundra above 3500m
X <- as.data.frame(scale(hta[,c(3:8,29)]))
y <- hta$Texpand
Xy <- as.data.frame(cbind(X,y))
bestHTE <- bestglm(Xy, IC = "AIC", family = binomial)
bestHTE # CV, Solar, Slope
bestHTE$BestModels
bestHTE <- glm(Texpand ~ CV + Solar + Slope, family = binomial, data = hta)
summary(bestHTE)
logisticPseudoR2s(bestHTE)
Dsquared(bestHTE, adjust = TRUE)
bestHTECV<-cv.glm(data=hta,glmfit=bestHTE,K=10)
bestHTECV$delta
bestHTEAUC<-AUC(model=bestHTE)
bestHTEAUC$AUC
# Null AIC 148.26, AIC 138.42
min(hta$CV)
max(hta$CV)
min(hta$Solar)
max(hta$Solar)
min(hta$Slope)
max(hta$Slope)

plot(vc$Mean, vc$CV) # Less snow is more variable
cor <- cor.test(vc$Mean, vc$CV, method = "pearson")
cor



############################### Random Forests ############################################
of.rf <- randomForest(as.factor(Oexpand) ~ Tree_Cover + Elevation + Solar + Slope + Mean + CV + Trend, data = loa, family = binomial(logit), ntree = 5000, importance = TRUE)
of.rf
importance(of.rf)
varImpPlot(of.rf)

sh.rf <- randomForest(as.factor(Sexpand) ~ Shrub_Cover + Elevation + Solar + Trend + Slope + Mean + CV, data = lsa, family = binomial(logit), ntree = 5000, importance = TRUE)
sh.rf
importance(sh.rf)
varImpPlot(sh.rf)

tu.rf <- randomForest(as.factor(Texpand) ~ CV + Solar + Slope + Tundra_Cover + Elevation + Mean + Trend, data = hta, family = binomial(logit), ntree = 5000, importance = TRUE)
tu.rf
importance(tu.rf)
varImpPlot(tu.rf)



########################## Geographically Weighted Regression ##############################
# With package GWmodel
# Run model selection manually following the gwr.model.selection protocol
# Then compare to the original bestglm model
loa <- SpatialPointsDataFrame(cbind(loa$X,loa$Y), loa)
DM <- gw.dist(dp.locat=cbind(loa$X,loa$Y))
bw.f2 <- bw.ggwr(Oexpand~Tree_Cover+Elevation+Solar,data=loa,dMat=DM,
                 family ="binomial")
res.binomial <- ggwr.basic(Oexpand~Tree_Cover+Elevation+Solar,bw=bw.f2,
                           data=loa,dMat=DM,family ="binomial")
res.binomial$GW.diagnostic

lsa <- SpatialPointsDataFrame(cbind(lsa$X,lsa$Y), lsa)
DM <- gw.dist(dp.locat=cbind(lsa$X,lsa$Y))
bw.f2 <- bw.ggwr(Sexpand~Shrub_Cover+Elevation+Solar+Trend,data=lsa,dMat=DM,
                 family ="binomial")
res.binomial <- ggwr.basic(Sexpand~Shrub_Cover+Elevation+Solar+Trend,bw=bw.f2,
                           data=lsa,dMat=DM,family ="binomial")
res.binomial$GW.diagnostic

hta <- SpatialPointsDataFrame(cbind(hta$X,hta$Y), hta)
DM <- gw.dist(dp.locat=cbind(hta$X,hta$Y))
bw.f2 <- bw.ggwr(Texpand~CV+Solar+Slope,data=hta,dMat=DM,
                 family ="binomial")
res.binomial <- ggwr.basic(Texpand~CV+Solar+Slope,bw=bw.f2,
                           data=hta,dMat=DM,family ="binomial")
res.binomial$GW.diagnostic



################################### Climate Data ###########################################
# Analysis and Figure 2
d <- read.csv("Summer.csv")
m <- lm(d$Summer_mean ~ d$Year)
summary(m)
# 1972 to 2008
da <- d[20:56,]
m1 <- lm(da$Summer_mean ~ da$Year)
summary(m1)
time_series <- ts(da[,2:length(da)], start = 1972, end = 2008)
m2 <- tslm(time_series ~ trend)
summary(m2)
# Figure 2
ggplot(da, aes(x=Year,y=Summer_mean)) + 
  geom_point(size = 3) +
  scale_x_continuous(name = "Year", breaks = seq(1972,2008,4)) +
  scale_y_continuous(name = "Mean Summer Temp. (˚C)", breaks = 
                       seq(3,12,1), 
                     limits = c(3,11)) + 
  geom_smooth(method=loess, color = "blue", fill = "blue", alpha = 0.1) +
  geom_smooth(method=lm, color = "red", "fill" = "red", alpha = 0.1) +
  theme_bw() +
  theme(legend.position="NULL",
        axis.title.x = element_text(face="bold", size = 18), 
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold",size=18))



############################ Bare to Tundra, Soil ##########################################
bt <- read.csv("BareTundra.csv", header = TRUE)
log <- bt[1:24,]
log$Change <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1)

BTlog<-glm(Change ~ Depth, family = binomial, data = log)
summary(BTlog) # NSD

BTlog<-glm(Change ~ pH, family = binomial, data = log)
summary(BTlog) # NSD

BTlog<-glm(Change ~ Bulk, family = binomial, data = log)
summary(BTlog) # NSD

########################### Tundra to Shrub, Soil ##########################################
ts <- read.csv("TundraShrub.csv")
log <- as.data.frame(ts[1:33,])
log$Change <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

TSlog<-glm(Change ~ Depth, family = binomial, data = log)
summary(TSlog) # NSD

TSlog<-glm(Change ~ pH, family = binomial, data = log)
summary(TSlog) # NSD

TSlog<-glm(Change ~ Bulk, family = binomial, data = log)
summary(TSlog) # NSD

############################# Tundra to Forest, Soil #######################################
log <- read.csv("TundraForest.csv")

# All Species
TOlog<-glm(Change ~ Depth, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Change ~ pH, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Change ~ Bulk, family = binomial, data = log)
summary(TOlog) # NSD

# Limber Pine
TOlog<-glm(Pine ~ Depth, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Pine ~ pH, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Pine ~ Bulk, family = binomial, data = log)
summary(TOlog) # NSD

# Englemann Spruce
TOlog<-glm(Spruce ~ Depth, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Spruce ~ pH, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Spruce ~ Bulk, family = binomial, data = log)
summary(TOlog) # Significant

# Subalpine Fir
TOlog<-glm(Fir ~ Depth, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Fir ~ pH, family = binomial, data = log)
summary(TOlog) # NSD

TOlog<-glm(Fir ~ Bulk, family = binomial, data = log)
summary(TOlog) # NSD
