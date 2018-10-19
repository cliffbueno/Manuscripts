# Analysis of fungal phenology from Niwot Ridge Saddle, Summer 2017
# by Cliff Bueno de Mesquita and Cormac Martinez del Rio
# Collected ArtSco, GenAlg, GeuRos, and PolBis throughout summer
# Looked at AMF, DSE, and FRE root colonization
# Make graphs, conduct repeated measures anova, test for effects of environment or plant variables
# Now published in Mycorrhiza

#################################### Setup #################################################
library(ggplot2)
library(plyr)
library(PMCMR)
library(car)
library(leaps)
library(bestglm)
library(AICcmodavg)
library(corrplot)
library(nlme)
library(rcompanion)
library(lme4)
library(lmerTest)
library(vegan)
source("~/Desktop/Functions/Summary.R")
setwd("~/Desktop/CU/2Research/FungalPhenology")
d <- read.csv("FungPhenAll.csv")
d$subject <- as.factor(d$subject)
mt <- read.csv("moist_temp.csv")
mt$subject <- as.factor(mt$subject)
mt$location <- as.factor(mt$location)
mt$temp2 <- mt$temp*2
mean(d$amf.hyphae)
mean(d$total.dse)
mean(d$fre.hyphae)
artsco <- subset(d, spp == "ArtSco")
genalg <- subset(d, spp == "GenAlg")
geuros <- subset(d, spp == "GeuRos")
polbis <- subset(d, spp == "PolBis")

############################ Moisture and Temperature ######################################
sumM <- summarySE(mt, measurevar ="moisture",groupvars=c("snowfree"))
sumT <- summarySE(mt, measurevar ="temp",groupvars=c("snowfree"))
sumT$temp2 <- sumT$temp*2
sumM$temp <- sumT$temp2
sumM$Tse <- sumT$se

# Figure 4
ggplot(sumM, aes(snowfree, moisture)) +
  geom_point(data=sumM,aes(snowfree, moisture), col = "blue", pch = 17) +
  geom_point(data=sumT,aes(snowfree, temp2), col = "red" ) +
  geom_errorbar(data=sumM,aes(ymin=moisture-se,ymax=moisture+se),width=.1,col= "blue") +
  geom_errorbar(data=sumM,aes(ymin=temp-Tse, ymax=temp+Tse), width=.1, col = "red") +
  geom_smooth(data=mt, aes(snowfree,moisture),method=loess,col="blue",se=FALSE,linetype=3) +
  geom_smooth(data=mt, aes(snowfree,temp2),method=loess,col="red",se=FALSE) +
  scale_y_continuous(sec.axis = sec_axis(~./2, name = "Soil Temperature (˚C)")) +
  xlab("Snow-free Days") +
  ylab("Soil Moisture (% VWC)") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size = 16, vjust = 0), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold", size=16))

# Biomass
sum <- summarySE(d, measurevar ="biomass",groupvars=c("snowfree","spp"))
pd <- position_dodge(0.2)
ggplot(sum, aes(snowfree, biomass)) +
  geom_errorbar(aes(ymin=biomass-se, ymax=biomass+se, colour = spp),width=.1, position=pd) +
  geom_line(position = pd, aes(shape = spp, linetype = spp, colour = spp, group = spp)) +
  geom_point(position = pd, size = 3, aes(shape = spp, linetype = spp, colour = spp, group = spp)) +
  scale_color_discrete("") +
  scale_shape_manual("",values = c(0,1,2,5)) +
  scale_linetype_manual("",values = c("solid","dashed","dotted","dotdash")) +
  xlab("Snow-free Days") +
  ylab("Biomass (g)") +
  labs(color = "Species") +
  theme_bw() +
  theme(legend.text = element_text(size = 18),
        legend.key.size = unit(3,"line"),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_blank())

model.a = gls(biomass ~ spp + snowfree + spp*snowfree,data=d)
ACF(model.a, form = ~ snowfree | subject)
model <- gls(biomass ~ spp + snowfree + spp*snowfree, correlation = corAR1(form = ~ snowfree | subject),data=d,method="REML")
Anova(model) # Significant Species and Snowfree
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

# Combined biomass/phenophase and structure graph, vertical stack, connected points
# Figure 3
g <-read.csv("structurepanel2.csv")
g$Structure = factor(g$Structure, levels=c("Dark.Septate","Coarse.Hyphae","Fine.Hyphae","Arbuscules","Vesicles"))
sum <- summarySE(g, measurevar = "Colonization",groupvars=c("snowfree","spp","Structure"))
pd <- position_dodge(1)
ggplot(sum, aes(snowfree, Colonization)) +
  geom_line(position = pd, aes(shape = spp, linetype = spp, colour = spp, group = spp)) +
  geom_point(position = pd, size = 3, aes(shape=spp, linetype=spp, colour = spp, group = spp)) +
  geom_errorbar(aes(ymin=Colonization-se,ymax=Colonization+se,colour=spp),width=.1,position=pd) +
  xlab("Snow-free Days") +
  ylab("% Colonization") +
  scale_color_discrete("") +
  scale_shape_manual("",values = c(0,1,2,5)) +
  scale_linetype_manual("",values = c("solid","dashed","dotted","dotdash")) +
  labs(color = "Species") +
  facet_wrap(~Structure, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_blank(),
        legend.background = element_rect(colour = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none" )



######################## Correlations between AMF, DSE, FRE ################################
plot(artsco$amf.hyphae, artsco$fre.hyphae)
plot(artsco$amf.hyphae, artsco$dse.hyphae)
plot(artsco$dse.hyphae, artsco$fre.hyphae)
cor.test(artsco$amf.hyphae, artsco$fre.hyphae) # Sig pos
cor.test(artsco$amf.hyphae, artsco$dse.hyphae) # Not sig
cor.test(artsco$dse.hyphae, artsco$fre.hyphae) # Not sig

plot(genalg$amf.hyphae, genalg$fre.hyphae)
plot(genalg$amf.hyphae, genalg$dse.hyphae)
plot(genalg$dse.hyphae, genalg$fre.hyphae)
cor.test(genalg$amf.hyphae, genalg$fre.hyphae) # Sig pos
cor.test(genalg$amf.hyphae, genalg$dse.hyphae) # Not sig
cor.test(genalg$dse.hyphae, genalg$fre.hyphae) # Not sig

plot(geuros$amf.hyphae, geuros$fre.hyphae)
plot(geuros$amf.hyphae, geuros$dse.hyphae)
plot(geuros$dse.hyphae, geuros$fre.hyphae)
cor.test(geuros$amf.hyphae, geuros$fre.hyphae) # Sig pos
cor.test(geuros$amf.hyphae, geuros$dse.hyphae) # Not sig
cor.test(geuros$dse.hyphae, geuros$fre.hyphae) # Not sig

plot(polbis$amf.hyphae, polbis$fre.hyphae)
plot(polbis$amf.hyphae, polbis$dse.hyphae)
plot(polbis$dse.hyphae, polbis$fre.hyphae)
cor.test(polbis$amf.hyphae, polbis$fre.hyphae) # Sig pos
cor.test(polbis$amf.hyphae, polbis$dse.hyphae) # Not sig
cor.test(polbis$dse.hyphae, polbis$fre.hyphae) # Not sig

line1 <- lm(artsco$fre.hyphae ~ artsco$amf.hyphae)
line2 <- lm(genalg$fre.hyphae ~ genalg$amf.hyphae)
line3 <- lm(geuros$fre.hyphae ~ geuros$amf.hyphae)
line4 <- lm(polbis$fre.hyphae ~ polbis$amf.hyphae)

# Figure 5
par(mfrow=c(4,3), oma=c(3,3,1,1) +0.1,mar=c(0.75,0.75,0.75,0.75) +0.1,mgp=c(2,1,0),xpd=NA)
par(xaxs="i")
plot(artsco$amf.hyphae, artsco$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); title(main="AMF vs. FRE", font.main=2, adj=0.5, cex.main=1.2, line=0.5); axis(side=2, at=c(0,20,40,60,80,100), labels=TRUE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); lines(artsco$amf.hyphae, fitted(line1), col="purple"); text(x=20, y=96, "rho=0.55, p<0.01", cex=0.75)
plot(artsco$amf.hyphae, artsco$dse.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); title(main="AMF vs. DSE", font.main=2, adj=0.5, cex.main=1.2, line=0.5); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.28, p>0.05", cex=0.75)
plot(artsco$dse.hyphae, artsco$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); title(main="DSE vs. FRE", font.main=2, adj=0.5, cex.main=1.2, line=0.5); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.20, p>0.05", cex=0.75)
plot(genalg$amf.hyphae, genalg$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=TRUE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); lines(genalg$amf.hyphae, fitted(line2), col="purple"); text(x=20, y=96, "rho=0.57, p<0.01", cex=0.75)
plot(genalg$amf.hyphae, genalg$dse.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=-0.17, p>0.05", cex=0.75)
plot(genalg$dse.hyphae, genalg$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.01, p>0.05", cex=0.75)
plot(geuros$amf.hyphae, geuros$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=TRUE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); lines(geuros$amf.hyphae, fitted(line3), col="purple"); text(x=20, y=96, "rho=0.52, p<0.01", cex=0.75)
plot(geuros$amf.hyphae, geuros$dse.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.17, p>0.05", cex=0.75)
plot(geuros$dse.hyphae, geuros$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n',xaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); axis(side=1, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.30, p>0.05", cex=0.75)
plot(polbis$amf.hyphae, polbis$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=TRUE); lines(polbis$amf.hyphae, fitted(line4), col="purple"); text(x=20, y=96, "rho=0.58, p<0.01", cex=0.75)
plot(polbis$amf.hyphae, polbis$dse.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.15, p>0.05", cex=0.75)
plot(polbis$dse.hyphae, polbis$fre.hyphae,ylab="",xlab="", cex.main=1, ylim=c(0,100),xlim=c(0,100),yaxt='n'); axis(side=2, at=c(0,20,40,60,80,100), labels=FALSE); text(x=20, y=96, "rho=0.05, p>0.05", cex=0.75)
mtext("Percent Root Colonization",side=2,outer=TRUE,cex=0.8,line=1.3)
mtext("Percent Root Colonization",side=1,outer=TRUE,cex=0.8,line=1.25)



############################### Variance Partitioning ######################################
# partition variance by environmental or plant variables for each species
# Results for Table 2
par(mfrow=c(1,1))
# Artsco
varp <- varpart(artsco$amf.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = artsco)
varp
plot (varp, digits = 2)
varp <- varpart(artsco$fre.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = artsco)
varp
plot (varp, digits = 2)
varp <- varpart(artsco$vesicle, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = artsco)
varp
plot (varp, digits = 2)
varp <- varpart(artsco$arbuscule, ~ phenophase + biomass, ~ mean.moist + mean.temp, data =artsco)
varp
plot (varp, digits = 2)
varp <- varpart(artsco$dse.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = artsco)
varp
plot (varp, digits = 2)

# Genalg
varp <- varpart(genalg$amf.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = genalg)
varp
plot (varp, digits = 2)
varp <- varpart(genalg$arbuscule, ~ phenophase + biomass, ~ mean.moist + mean.temp, data =genalg)
varp
plot (varp, digits = 2)
varp <- varpart(genalg$vesicle, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = genalg)
varp
plot (varp, digits = 2)
varp <- varpart(genalg$fre.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = genalg)
varp
plot (varp, digits = 2)
varp <- varpart(genalg$dse.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = genalg)
varp
plot (varp, digits = 2)

# Geuros
varp <- varpart(geuros$amf.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = geuros)
varp
plot (varp, digits = 2)
varp <- varpart(geuros$arbuscule, ~ phenophase + biomass, ~ mean.moist + mean.temp, data =geuros)
varp
plot (varp, digits = 2)
varp <- varpart(geuros$vesicle, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = geuros)
varp
plot (varp, digits = 2)
varp <- varpart(geuros$fre.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = geuros)
varp
plot (varp, digits = 2)
varp <- varpart(geuros$dse.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = geuros)
varp
plot (varp, digits = 2)

# Polbis
varp <- varpart(polbis$amf.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = polbis)
varp
plot (varp, digits = 2)
varp <- varpart(polbis$arbuscule, ~ phenophase + biomass, ~ mean.moist + mean.temp, data =polbis)
varp
plot (varp, digits = 2)
varp <- varpart(polbis$vesicle, ~ phenophase + biomass, ~ mean.moist + mean.temp, data = polbis)
varp
plot (varp, digits = 2)
varp <- varpart(polbis$fre.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp, data=polbis)
varp
plot (varp, digits = 2)
varp <- varpart(polbis$dse.hyphae, ~ phenophase + biomass, ~ mean.moist + mean.temp,data= polbis)
varp
plot (varp, digits = 2)


################### Repeated Measures ANOVA, Time*Species, all data ########################
# Do for each fungal structure
# Results shown on Figure 3
model <- lme(amf.hyphae ~ spp + snowfree + spp*snowfree, random = ~1|observer, correlation = corAR1(form = ~ snowfree | observer/subject),data=d,method="REML")
Anova(model) # Significant Species and Snowfree
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model <- lme(fre.hyphae ~ spp + snowfree + spp*snowfree, random = ~1|observer, correlation = corAR1(form = ~ snowfree | observer/subject),data=d,method="REML")
Anova(model) # Significant Species, Snowfree and Interaction
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model <- lme(arbuscule ~ spp + snowfree + spp*snowfree, random = ~1|observer, correlation = corAR1(form = ~ snowfree | observer/subject),data=d,method="REML")
Anova(model) # Significant Species, Snowfree
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model <- lme(vesicle ~ spp + snowfree + spp*snowfree, random = ~1|observer, correlation = corAR1(form = ~ snowfree | observer/subject),data=d,method="REML")
Anova(model) # Significant Species, Snowfree and Interaction
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model <- lme(dse.hyphae ~ spp + snowfree + spp*snowfree, random = ~1|observer, correlation = corAR1(form = ~ snowfree | observer/subject),data=d,method="REML")
Anova(model) # NS
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))



####################### Multivariate Models within Species #################################
# phenophase, biomass, moist, temp
# Remove unused phenophase levels for artsco and genalg
# Results for Table 1

### ArtSco
toremove <- which(d$phenophase=="budding"|d$spp=="genalg"|d$spp=="geuros"|d$spp=="polbis")
artscoDL <- droplevels(d[-toremove,])
model<-lmer(amf.hyphae~phenophase+biomass+mean.moist+mean.temp+(1|observer),data=artscoDL)
Anova(model) # Significant biomass, moist and temp
model<-lmer(fre.hyphae~phenophase+biomass+mean.moist+mean.temp+(1|observer),data=artscoDL)
Anova(model) # Significant temp and biomass
model <- lmer(arbuscule~phenophase+biomass+mean.moist+mean.temp+(1|observer), data=artscoDL)
Anova(model) # Significant moist
model <- lmer(vesicle ~ phenophase+biomass+mean.moist+mean.temp+(1|observer), data=artscoDL)
Anova(model) # Nothing significant
model<-lmer(dse.hyphae~phenophase+biomass+mean.moist+ mean.temp+(1|observer), data=artscoDL)
Anova(model) # Nothing significant

### GenAlg
toremove<-which(d$phenophase=="flowering"|d$spp=="artsco"|d$spp=="geuros"|d$spp=="polbis")
genalgDL <- droplevels(d[-toremove,])
model<-lmer(amf.hyphae~phenophase+biomass+mean.moist+mean.temp+(1|observer),data=genalgDL)
Anova(model) # Significant biomass, marginal temp
model<-lmer(fre.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer),data=genalgDL)
Anova(model) # Significant phenophase
model<-lmer(arbuscule~phenophase+biomass+mean.moist+mean.temp+(1|observer), data=genalgDL)
Anova(model) # Phenophase, moist, temp significant
model <- lmer(vesicle~phenophase+biomass+mean.moist+mean.temp+(1|observer), data=genalgDL)
Anova(model) # Nothing significant
model<-lmer(dse.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer),data=genalgDL)
Anova(model) # Nothing significant 

### GeuRos
model <-lmer(amf.hyphae~phenophase+biomass+mean.moist+mean.temp+(1|observer), data=geuros)
Anova(model) # All significant
model <-lmer(fre.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer), data=geuros)
Anova(model) # Phenophase, temp, moist significant
model <- lmer(arbuscule ~ phenophase+biomass+mean.moist+mean.temp+(1|observer), data=geuros)
Anova(model) # Phenophase significant, temp, moist marginal
model <- lmer(vesicle ~ phenophase+biomass+mean.moist+mean.temp+(1|observer), data=geuros)
Anova(model) # Nothing significant
model <- lmer(dse.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer),data=geuros)
Anova(model) # Phenophase, temp, moist significant

### PolBis
model<-lmer(amf.hyphae~phenophase+biomass+mean.moist+mean.temp+(1|observer),data=polbis)
Anova(model) #  Phenophase, moist, temp significant
model<-lmer(fre.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer),data=polbis)
Anova(model) # Nothing significant
model <- lmer(arbuscule ~ phenophase+biomass+mean.moist+mean.temp+(1|observer), data=polbis)
Anova(model) # Nothing significant
model <- lmer(vesicle ~ phenophase+biomass+mean.moist+mean.temp + (1|observer), data=polbis)
Anova(model) # Nothing significant
model <- lmer(dse.hyphae ~ phenophase+biomass+mean.moist+mean.temp+(1|observer),data=polbis)
Anova(model) # Biomass significant
