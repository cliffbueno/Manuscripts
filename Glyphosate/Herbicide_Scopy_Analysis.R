# Analysis of fungal colonization in plant roots from Chatfield Farm
# Plants were sampled in June 2019 by CB, AS, CM and AT
# Analysis by Cliff Bueno de Mesquita, September 2019

################################# Setup ############################################
library(ggplot2)
library(Matrix)
library(lme4)
library(PMCMR)
library(lmodel2)
library(car)
library(nlme)
library(rcompanion)
library(emmeans)
library(pscl)
library(reshape2)
library(multcomp)
library(FSA)
source("~/Desktop/OneDrive - UCB-O365/Functions/Summary.R")
setwd("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide")
data <- read.csv("Herbicide_Scopy_Data.csv")
data$Herbicide <- as.factor(data$Herbicide)
table(data$Species)
data <- subset(data, Code == "BrIn" | Code == "EuCh")
b <- subset(data, Code == "BrIn")
table(b$Herbicide)
e <- subset(data, Code == "EuCh")
table(e$Herbicide)
min(b$Tot.Perc.Col)
max(b$Tot.Perc.Col)
mean(b$Tot.Perc.Col)
se(b$Tot.Perc.Col)
mean(e$Tot.Perc.Col)
se(e$Tot.Perc.Col)

################################## Analysis ########################################
model <- lmer(Tot.Perc.Col ~ Herbicide + (1|Viewer), data = b)
Anova(model) # Sig
emmeans(model, list(pairwise ~ Herbicide), adjust = "tukey")
m <- aov(Tot.Perc.Col ~ Herbicide, data = b)
summary(m)
TukeyHSD(m)
model <- lmer(Tot.Perc.Col ~ Herbicide + (1|Viewer), data = e)
Anova(model) # NS

ggplot(data,aes(x = Herbicide, y = Tot.Perc.Col, fill = Species)) + 
  geom_boxplot() + 
  labs(y = "% Root Colonization",
       x = "# Glyphosate Applications") +
  theme_bw() +
  theme(legend.position = c(0.77,0.8),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

ggplot(b,aes(x = Herbicide, y = Tot.Perc.Col)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(size = 4, alpha = 0.5) +
  geom_text(aes(x = 1, y = 45, label = "a"), size = 5, colour = "black") +
  geom_text(aes(x = 2, y = 45, label = "ab"), size = 5, colour = "black") +
  geom_text(aes(x = 3, y = 45, label = "b"), size = 5, colour = "black") +
  labs(y = "% Root Colonization",
       x = "# Herbicide Applications") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

model <- lmer(FRE.Tot ~ Herbicide + (1|Viewer), data = b)
Anova(model) # NS
model <- lmer(FRE.Tot ~ Herbicide + (1|Viewer), data = e)
Anova(model) # NS

ggplot(data,aes(x = Herbicide, y = FRE.Tot, fill = Species)) + 
  geom_boxplot() + 
  labs(y = "% FRE Colonization",
       x = "# Glyphosate Applications") +
  theme_bw() +
  theme(legend.position = c(0.77,0.8),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

model <- lmer(AMF.Tot ~ Herbicide + (1|Viewer), data = b)
Anova(model) # Marginal
emmeans(model, list(pairwise ~ Herbicide), adjust = "tukey")
m <- aov(AMF.Tot ~ Herbicide, data = b)
summary(m)
TukeyHSD(m)
model <- lmer(AMF.Tot ~ Herbicide + (1|Viewer), data = e)
Anova(model) # NS

ggplot(data,aes(x = Herbicide, y = AMF.Tot, fill = Species)) + 
  geom_boxplot() + 
  labs(y = "% AMF Colonization",
       x = "# Glyphosate Applications") +
  theme_bw() +
  theme(legend.position = c(0.77,0.8),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

model <- lmer(DSE.Tot ~ Herbicide + (1|Viewer), data = b)
Anova(model) # NS
model <- lmer(DSE.Tot ~ Herbicide + (1|Viewer), data = e)
Anova(model) # NS

ggplot(data,aes(x = Herbicide, y = DSE.Tot, fill = Species)) + 
  geom_boxplot() + 
  labs(y = "% DSE Colonization",
       x = "# Glyphosate Applications") +
  theme_bw() +
  theme(legend.position = c(0.77,0.8),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

# Make huge facet grid with plant by fungus panels
# Expand data
long <- melt(data, id.vars = c("Code","Herbicide"), measure.vars = c("Tot.Perc.Col","AMF.Tot","DSE.Tot","FRE.Tot"))
spec_names <- c('BrIn'="Bromus inermis",'EuCh'="Eucrypta chrysanthemifolia")
fung_names <- c('Tot.Perc.Col'="Total",'AMF.Tot'="AMF",'DSE.Tot'="DSE",'FRE.Tot'="FRE")
label.df <- data.frame(
  Code = c("BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn","BrIn",
           "EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh","EuCh"),
  variable = c("Tot.Perc.Col","Tot.Perc.Col","Tot.Perc.Col","AMF.Tot","AMF.Tot","AMF.Tot",
               "DSE.Tot","DSE.Tot","DSE.Tot","FRE.Tot","FRE.Tot","FRE.Tot",
               "Tot.Perc.Col","Tot.Perc.Col","Tot.Perc.Col","AMF.Tot","AMF.Tot","AMF.Tot",
               "DSE.Tot","DSE.Tot","DSE.Tot","FRE.Tot","FRE.Tot","FRE.Tot"),
  Herbicide = c("0","2","4","0","2","4","0","2","4","0","2","4",
                "0","2","4","0","2","4","0","2","4","0","2","4"),
  Value = c(45,45,45,45,45,45,45,45,45,45,45,45,
            45,45,45,45,45,45,45,45,45,45,45,45),
  Sig = c("a","ab","b","a","ab","b"," ","NSD"," "," ","NSD"," ",
          " ","NSD"," "," ","NSD"," "," ","NSD"," "," ","NSD"," "))
ggplot(long,aes(x = Herbicide, y = value)) + 
  geom_boxplot() +
  geom_text(data = label.df, aes(x = Herbicide,y=Value,label=Sig,group=NULL)) +
  labs(y = "% Root Colonization",
       x = "# Glyphosate Applications") +
  facet_grid(variable ~ Code, labeller = labeller(Code = spec_names,variable = fung_names)) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18),
        strip.text.x = element_text(size = 14, face = "italic"),
        strip.text.y = element_text(size = 14)) +
  ylim(0,47)



############################  Binomial Models ######################################
m <- glm(cbind(AMF.Yes, AMF.No) ~ 0 + Herbicide, family = binomial, data = b)
summary(m)
Anova(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey"))) # a, b, c
plogis(coef(m))
plogis(confint(m))

m <- glm(cbind(FRE.Yes, FRE.No) ~ 0 + Herbicide, family = binomial, data = b)
summary(m)
Anova(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey"))) # NSD
plogis(coef(m))
plogis(confint(m))

m <- glm(cbind(DSE.Yes, DSE.No) ~ 0 + Herbicide, family = binomial, data = b)
summary(m)
Anova(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey"))) # a, a, b
plogis(coef(m))
plogis(confint(m))

m <- glm(cbind(AMF.Yes, AMF.No) ~ 0 + Herbicide, family = binomial, data = e)
summary(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey")))
plogis(coef(m))
plogis(confint(m))

m <- glm(cbind(FRE.Yes, FRE.No) ~ 0 + Herbicide, family = binomial, data = e)
summary(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey")))
plogis(coef(m))
plogis(confint(m))

m <- glm(cbind(DSE.Yes, DSE.No) ~ 0 + Herbicide, family = binomial, data = e)
summary(m)
summary(glht(m,linfct=mcp(Herbicide="Tukey")))
plogis(coef(m))
plogis(confint(m))

# Graph for paper
long <- melt(data, id.vars = c("Code","Herbicide"), measure.vars = c("Tot.Perc.Col","AMF.Tot","DSE.Tot","FRE.Tot"))
spec_names <- c('BrIn'="Bromus inermis",'EuCh'="Eucrypta chrysanthemifolia")
fung_names <- c('Tot.Perc.Col'="Total",'AMF.Tot'="a) AMF",'DSE.Tot'="b) DSE",'FRE.Tot'="c) FRE")
long <- subset(long, Code == "BrIn")
long <- subset(long, variable != "Tot.Perc.Col")
ci <- read.csv("scopy_confints.csv")
ci$Herbicide <- as.factor(ci$Herbicide)
ci$label <- as.character(ci$label)
ggplot(long, aes(x=Herbicide, y=value)) + 
  geom_jitter(size = 2, alpha = 0.3,
              position=position_jitter(width=0.1, height=0.1)) +
  geom_segment(aes(xend=Herbicide, y=lo, yend=hi),
               data=ci, size=3, col="blue", alpha=.3) + 
  geom_point(aes(x=Herbicide, y=prob), data=ci, col="blue", size=2) +
  geom_text(data = ci, aes(x = Herbicide, y = value, label = label), size = 5) +
  labs(y = "% Root Colonization",
       x = "# Herbicide Applications") +
  facet_wrap(~ variable, labeller = as_labeller(fung_names)) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18),
        strip.text = element_text(size = 14))
