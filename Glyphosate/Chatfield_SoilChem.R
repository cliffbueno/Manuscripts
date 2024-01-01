# Chatfield Soil Chemistry Analysis

#### Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(Matrix)
library(lme4)
library(PMCMRplus)
library(lmodel2)
library(car)
library(nlme)
library(rcompanion)
library(emmeans)
library(pscl)
library(reshape2)
library(cowplot)
library(vegan)
library(RVAideMemoire)
library(phyloseq)
library(multcomp)

# Plot theme
theme_set(theme_bw())

# Working directory
setwd("~/Desktop/OneDrive - UCB-O365/CU/2Research/Herbicide")

# Data
d <- read.csv("Chatfield_SoilChem.csv")
d$Time_point <- as.factor(d$Time_point)
d$plot <- as.factor(d$plot)
table(d$Time_point)
table(d$Treatment)
d$soil_type
mean(d$texture_sand_percent, na.rm = T)
mean(d$texture_silt_percent, na.rm = T)
mean(d$texture_clay_percent, na.rm = T)
mean(d$organic_matter_percent)
c <- subset(d, Treatment == "0H")
mean(c$pH)

# Update - try column with number of applications actually applied before the sampling
d$Applications <- NA
for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "0H") {
    d$Applications[i] <- 0
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "2H" & d$Time_point[i] == 1) {
    d$Applications[i] <- 1
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "2H" & d$Time_point[i] == 2) {
    d$Applications[i] <- 2
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "2H" & d$Time_point[i] == 3) {
    d$Applications[i] <- 2
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "2H" & d$Time_point[i] == 4) {
    d$Applications[i] <- 2
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "5H" & d$Time_point[i] == 1) {
    d$Applications[i] <- 1
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "5H" & d$Time_point[i] == 2) {
    d$Applications[i] <- 2
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "5H" & d$Time_point[i] == 3) {
    d$Applications[i] <- 4
  }
}

for (i in 1:nrow(d)) {
  if (d$Treatment[i] == "5H" & d$Time_point[i] == 4) {
    d$Applications[i] <- 5
  }
}
d$Applications <- as.factor(d$Applications)

table(d$Applications)
# Analyses - For each variable, test assumptions and effects of time and herbicide
# Variables - OM, P weak, P strong, K, Mg, Ca, pH, CEC, Nitrate, Ammonium, Total N

# Contrasts are important for Type III ANOVA
# Default is treatment. Define Helmert here
my.contrasts <- list(Treatment = "contr.Helmert", Time_point = "contr.Helmert")

# OM
leveneTest(d$organic_matter_percent ~ d$Treatment)
leveneTest(d$organic_matter_percent ~ d$Time_point)
m <- aov(d$organic_matter_percent ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time
Anova(m, type = "II") # NSD
Anova(m, type = "III") # NSD
m <- lmer(organic_matter_percent ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # NSD
Anova(m, type = "II") # Time
m <- lmer(organic_matter_percent ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(organic_matter_percent ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II") # NSD

m <- lmer(organic_matter_percent ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# P weak (currently available to plants)
leveneTest(d$P1_weakBray_ppm ~ d$Treatment)
leveneTest(d$P1_weakBray_ppm ~ d$Time_point)
m <- aov(d$P1_weakBray_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time
Anova(m, type = "II") # Time
Anova(m, type = "III") # Time
m <- lmer(P1_weakBray_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time
Anova(m, type = "II") # Time
m <- lmer(P1_weakBray_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(P1_weakBray_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II") # Time

m <- lmer(P1_weakBray_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# P strong (reserve of P)
leveneTest(d$P2_strongBray_ppm ~ d$Treatment)
leveneTest(d$P2_strongBray_ppm ~ d$Time_point)
m <- aov(d$P2_strongBray_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time
Anova(m, type = "II") # NSD
Anova(m, type = "III") # NSD
m <- lmer(P2_strongBray_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time marginal
Anova(m, type = "II") # Time
m <- lmer(P2_strongBray_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(P2_strongBray_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II") # Time

m <- lmer(P2_strongBray_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# K
leveneTest(d$K_ppm ~ d$Treatment)
leveneTest(d$K_ppm ~ d$Time_point)
m <- aov(d$K_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # NSD
Anova(m, type = "II") # NSD
Anova(m, type = "III") # NSD
m <- lmer(K_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # NSD
Anova(m, type = "II") # Time
m <- lmer(K_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(K_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")

m <- lmer(K_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# Mg
leveneTest(d$Mg_ppm ~ d$Treatment)
leveneTest(d$Mg_ppm ~ d$Time_point)
m <- aov(d$Mg_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # NSD
Anova(m, type = "II") # NSD
Anova(m, type = "III") # NSD
m <- lmer(Mg_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time
Anova(m, type = "II") # Time
m <- lmer(Mg_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(Mg_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")

m <- lmer(Mg_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# Ca
leveneTest(d$Ca_ppm ~ d$Treatment)
leveneTest(d$Ca_ppm ~ d$Time_point)
m <- aov(d$Ca_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time and Herbicide
TukeyHSD(m) # 5 different than 0 and 2!
Anova(m, type = "II") # Time and Herbicide
Anova(m, type = "III") # Herbicide
m <- lmer(Ca_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time marginal, Herbicide
Anova(m, type = "II") # Time, Herbicide
m <- lmer(Ca_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
summary(glht(m,linfct=mcp(Treatment="Tukey")))
m <- lmer(Ca_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")
summary(glht(m,linfct=mcp(Treatment="Tukey")))

m <- lmer(Ca_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# pH
leveneTest(d$pH ~ d$Treatment)
leveneTest(d$pH ~ d$Time_point)
m <- aov(d$pH ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Herbicide
TukeyHSD(m) # All different
Anova(m, type = "II") # Herbicide
Anova(m, type = "III") # Herbicide
m <- lmer(pH ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time marginal, Herbicide
Anova(m, type = "II") # Time
m <- lmer(pH ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
summary(glht(m,linfct=mcp(Treatment="Tukey")))
m <- lmer(pH ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")
summary(glht(m,linfct=mcp(Treatment="Tukey")))

m <- lmer(pH ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")
summary(glht(m,linfct=mcp(Applications="Tukey")))

# CEC
leveneTest(d$cation_exchange_meq_per_100g ~ d$Treatment)
leveneTest(d$cation_exchange_meq_per_100g ~ d$Time_point)
m <- aov(d$cation_exchange_meq_per_100g ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time and Herbicide
TukeyHSD(m) # 0 and 5
Anova(m, type = "II") # Same
Anova(m, type = "III") # Herbicide marginal
m <- lmer(cation_exchange_meq_per_100g ~ Treatment * Time_point + (1|plot),
          data = d)
summary(glht(m,linfct=mcp(Treatment="Tukey")))
Anova(m, type = "III") # Time, Herbicide, marginal interaction
Anova(m, type = "II") # Time
m <- lmer(cation_exchange_meq_per_100g ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
summary(glht(m,linfct=mcp(Treatment="Tukey")))
m <- lmer(cation_exchange_meq_per_100g ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II") # Time

m <- lmer(cation_exchange_meq_per_100g ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# Nitrate
leveneTest(d$Nitrate_ppm ~ d$Treatment)
leveneTest(d$Nitrate_ppm ~ d$Time_point)
m <- aov(d$Nitrate_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # Time and Herbicide and Interaction
TukeyHSD(m) # All different
Anova(m, type = "II") # Same
Anova(m, type = "III") # Time and Herbicide and Interaction
m <- lmer(Nitrate_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time and Herbicide and Interaction
Anova(m, type = "II") # Time and Herbicide and Interaction
summary(glht(m,linfct=mcp(Treatment="Tukey")))
m <- lmer(Nitrate_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time and Herbicide and Interaction
summary(glht(m,linfct=mcp(Treatment="Tukey")))
m <- lmer(Nitrate_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")

m <- lmer(Nitrate_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")
summary(glht(m,linfct=mcp(Applications="Tukey")))

# Ammonia
leveneTest(d$Ammonia_N_ppm ~ d$Treatment)
leveneTest(d$Ammonia_N_ppm ~ d$Time_point)
m <- aov(d$Ammonia_N_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # NSD
Anova(m, type = "II") # Same
Anova(m, type = "III") # NSD
m <- lmer(Ammonia_N_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # Time marginal
Anova(m, type = "II") # NSD
m <- lmer(Ammonia_N_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(Ammonia_N_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")

m <- lmer(Ammonia_N_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")

# Total N
leveneTest(d$Total_N_ppm ~ d$Treatment)
leveneTest(d$Total_N_ppm ~ d$Time_point)
m <- aov(d$Total_N_ppm ~ d$Treatment*d$Time_point)
shapiro.test(m$residuals)
summary(m) # NSD
Anova(m, type = "II") # Same
Anova(m, type = "III") # NSD
m <- lmer(Total_N_ppm ~ Treatment * Time_point + (1|plot),
          data = d)
Anova(m, type = "III") # NSD
Anova(m, type = "II") # NSD
m <- lmer(Total_N_ppm ~ Treatment * Time_point + (1|plot),
          contrasts = my.contrasts,
          data = d)
Anova(m, type = "III") # Time
m <- lmer(Total_N_ppm ~ Treatment + Time_point + (1|plot),
          data = d)
Anova(m, type = "II")

m <- lmer(Total_N_ppm ~  Date + Applications + (1|plot),
          data = d)
Anova(m, type = "II")



#### Final Timepoint ####
# The final sampling is the only true time when it's actually 0 vs. 2 vs. 5 sprays received...
t4 <- subset(d, Time_point == 4)

summary(aov(organic_matter_percent ~ Treatment, data = t4))
summary(aov(P1_weakBray_ppm ~ Treatment, data = t4))
summary(aov(P2_strongBray_ppm ~ Treatment, data = t4))
summary(aov(K_ppm ~ Treatment, data = t4))
summary(aov(Mg_ppm ~ Treatment, data = t4))
summary(aov(Ca_ppm ~ Treatment, data = t4))
summary(aov(pH ~ Treatment, data = t4))
TukeyHSD(aov(pH ~ Treatment, data = t4))
summary(aov(cation_exchange_meq_per_100g ~ Treatment, data = t4))
summary(aov(Nitrate_ppm ~ Treatment, data = t4))
TukeyHSD(aov(Nitrate_ppm ~ Treatment, data = t4))
summary(aov(Ammonia_N_ppm ~ Treatment, data = t4))



#### Significant ####
# Significant effects of herbicide on Ca, pH, CEC, NO3. So, make a 4 panel graph
# Ca
sumCa <- ddply(d, c("Treatment","Time_point"), summarise,
               meanCa <- mean(Ca_ppm),
               seCa = se(Ca_ppm))
colnames(sumCa) <- c("Treatment","Time_point","meanCa","seCa")
pd <- position_dodge(0.3)
g1<-ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = Ca_ppm), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = Ca_ppm, colour = Time_point),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumCa, aes(x = Treatment, ymin=meanCa-seCa, ymax=meanCa+seCa,
                                  colour = Time_point),width=.4, position=pd, size =0.9) +
  geom_point(data = sumCa, aes(Treatment, meanCa, colour = Time_point), 
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 3600, label = "Herbicide p = 0.015\nTime p < 0.001"), 
            size = 4, colour = "black") +
  geom_text(aes(x = 1, y = 3250, label = "a"), size = 4, colour = "black") +
  geom_text(aes(x = 2, y = 3250, label = "ab"), size = 4, colour = "black") +
  geom_text(aes(x = 3, y = 3250, label = "b"), size = 4, colour = "black") +
  labs(x = "Herbicide Applications", y = "[Calcium] (ppm)", title = "c) Calcium") +
  ylim(1750,3700) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

# pH
sumpH <- ddply(d, c("Treatment","Time_point"), summarise,
               meanpH <- mean(pH),
               sepH = se(pH))
colnames(sumpH) <- c("Treatment","Time_point","meanpH","sepH")
pd <- position_dodge(0.3)
g2<-ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = pH), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = pH, colour = Time_point),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumpH, aes(x = Treatment, ymin=meanpH-sepH, ymax=meanpH+sepH,
                                  colour = Time_point),width=.4, position=pd, size =0.9) +
  geom_point(data = sumpH, aes(Treatment, meanpH, colour = Time_point), 
             size = 4, position = pd) +
  geom_text(aes(x = 3, y = 7.55, label = "Herbicide p = 0.005\nTime p < 0.001"), 
            size = 4, colour = "black") +
  geom_text(aes(x = 1, y = 7.35, label = "a"), size = 4, colour = "black") +
  geom_text(aes(x = 2, y = 7.35, label = "ab"), size = 4, colour = "black") +
  geom_text(aes(x = 3, y = 7.35, label = "b"), size = 4, colour = "black") +
  labs(x = "Herbicide Applications", y = "pH", title = "b) pH") +
  ylim(6.4, 7.6) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

# CEC
sumCEC <- ddply(d, c("Treatment","Time_point"), summarise,
               meanCEC <- mean(cation_exchange_meq_per_100g),
               seCEC = se(cation_exchange_meq_per_100g))
colnames(sumCEC) <- c("Treatment","Time_point","meanCEC","seCEC")
pd <- position_dodge(0.3)
g3<-ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = cation_exchange_meq_per_100g),
               alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = cation_exchange_meq_per_100g, colour = Time_point),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumCEC, aes(x = Treatment, ymin=meanCEC-seCEC, ymax=meanCEC+seCEC,
                                  colour = Time_point),width=.4, position=pd, size =0.9) +
  geom_point(data = sumCEC, aes(Treatment, meanCEC, colour = Time_point), 
             size = 4, position = pd) +
  geom_text(aes(x = 2.75, y = 25.5, 
                label = "Herbicide p = 0.31\nTime p < 0.001\nHerbicide x Time p = 0.05"), 
            size = 4, colour = "black") +
#  geom_text(aes(x = 1, y = 23, label = "a"), size = 4, colour = "black") +
#  geom_text(aes(x = 2, y = 23, label = "a"), size = 4, colour = "black") +
#  geom_text(aes(x = 3, y = 23, label = "b"), size = 4, colour = "black") +
  labs(y = "CEC (meq/100g)", title = "d) Cation Exchange Capacity") +
  ylim(14,26.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

# NO3
sumNO3 <- ddply(d, c("Treatment","Time_point"), summarise,
                meanNO3 <- mean(Nitrate_ppm),
                seNO3 = se(Nitrate_ppm))
colnames(sumNO3) <- c("Treatment","Time_point","meanNO3","seNO3")
pd <- position_dodge(0.3)
g4<-ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = Nitrate_ppm),
               alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = Nitrate_ppm, colour = Time_point),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumNO3, aes(x = Treatment, ymin=meanNO3-seNO3, ymax=meanNO3+seNO3,
                                   colour = Time_point),width=.4, position=pd, size =0.9) +
  geom_point(data = sumNO3, aes(Treatment, meanNO3, colour = Time_point), 
             size = 4, position = pd) +
  geom_text(aes(x = 2.75, y = 37,
            label = "Herbicide p < 0.001\nTime p < 0.001\nHerbicide x Time p = 0.003"), 
            size = 4, colour = "black") +
  geom_text(aes(x = 1, y = 30, label = "a"), size = 4, colour = "black") +
  geom_text(aes(x = 2, y = 30, label = "b"), size = 4, colour = "black") +
  geom_text(aes(x = 3, y = 30, label = "c"), size = 4, colour = "black") +
  labs(y = "[Nitrate] (ppm)", title = "a) Nitrate") +
  ylim(0,41) +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 14), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

multiplot <- plot_grid(g4,g2,g1,g3, ncol=2, nrow=2, align = "hv")
multiplot # 13.61 x 7.06
pdf("Figure2_forPPT.pdf", width = 9, height = 6)
multiplot
dev.off()

# Make extra nitrate graph showing interaction
# Just June samplings and show line.
# NO3
j <- subset(d, Time_point == "1" | Time_point == "4")
sumNO3 <- ddply(j, c("Treatment","Time_point"), summarise,
                meanNO3 <- mean(Nitrate_ppm),
                seNO3 = se(Nitrate_ppm))
colnames(sumNO3) <- c("Treatment","Time_point","meanNO3","seNO3")
pd <- position_dodge(0.3)
ggplot() +
  geom_point(data = j,aes(x = Treatment, y = Nitrate_ppm, colour = Time_point),
             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumNO3, aes(x = Treatment, ymin=meanNO3-seNO3, ymax=meanNO3+seNO3,
                                   colour = Time_point),width=.4, position=pd, size =0.9) +
  geom_point(data = sumNO3, aes(Treatment, meanNO3, colour = Time_point), 
             size = 4, position = pd) +
  geom_line(data = sumNO3, aes(Treatment, meanNO3, colour = Time_point, group = Time_point), 
             size = 1, position = pd) +
  scale_colour_discrete(labels = c("June 2018","June 2019")) +
  scale_x_discrete(labels = c("0","2","5")) +
  labs(y = "[Nitrate] (ppm)", x = "# Herbicide Applications") +
  theme(legend.position = c(0.8,0.2),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16, vjust = -0.5),
        axis.title = element_text(face="bold", size = 16), 
        axis.text = element_text(size = 14),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
# Or facet wrap
j.long <- melt(j, id.vars = c("Treatment","Time_point"), measure.vars = "Nitrate_ppm")
facet_names <-  c(`1` = "a) June 2018",`4` = "b) June 2019")
ggplot() +
  geom_boxplot(data = j.long,aes(x = Treatment, y = value), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = j.long,aes(x = Treatment, y = value),
             size = 3, alpha = 0.2) +
  labs(x = "# Glyphosate Applications/Year", y = "[Nitrate] (ppm)") +
  scale_x_discrete(labels = c("0","2","5")) +
  facet_wrap(~ Time_point, labeller = as_labeller(facet_names)) +
  theme(axis.title = element_text(face="bold", size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))

# For Steve, for report, just boxplots of Calcium and pH
# Ca
g1 <- ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = Ca_ppm), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = Ca_ppm),
             size = 3, alpha = 0.2) +
  geom_text(aes(x = 3, y = 3600, label = "Herbicide p = 0.0001\nTime p = 0.001"), 
            size = 5, colour = "black") +
  geom_text(aes(x = 1, y = 3270, label = "a"), size = 5, colour = "black") +
  geom_text(aes(x = 2, y = 3270, label = "ab"), size = 5, colour = "black") +
  geom_text(aes(x = 3, y = 3270, label = "b"), size = 5, colour = "black") +
  labs(x = "Herbicide Applications", y = "[Calcium] (ppm)", title = "a) Calcium") +
  ylim(1750,3700) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 16), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
# pH
g2 <- ggplot() +
  geom_boxplot(data = d,aes(x = Treatment, y = pH), alpha = 0.5, outlier.shape = NA) +
  geom_point(data = d,aes(x = Treatment, y = pH),
             size = 3, alpha = 0.2) +
  geom_text(aes(x = 3, y = 7.55, label = "Herbicide p < 0.0001\nTime p = 0.11"), 
            size = 5, colour = "black") +
  geom_text(aes(x = 1, y = 7.37, label = "a"), size = 5, colour = "black") +
  geom_text(aes(x = 2, y = 7.37, label = "b"), size = 5, colour = "black") +
  geom_text(aes(x = 3, y = 7.37, label = "c"), size = 5, colour = "black") +
  labs(x = "Herbicide Applications", y = "pH", title = "b) pH") +
  ylim(6.4, 7.6) +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, vjust = -0.5),
        axis.title.y = element_text(face="bold", size = 16), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

multiplot <- plot_grid(g1,g2, ncol=1, align = "hv")
multiplot # 6.88 x 6.49



#### Realized Apps ####
d <- d %>%
  mutate(Date = recode_factor(Time_point,
                              "1" = "June 2018",
                              "2" = "August 2018",
                              "3" = "April 2019",
                              "4" = "June 2019"))
ggplot(d, aes(Date, Nitrate_ppm, colour = Applications)) +
  geom_point()

sumNO3 <- ddply(d, c("Date", "Treatment", "Applications"), summarise,
                meanNO3 <- mean(Nitrate_ppm),
                seNO3 = se(Nitrate_ppm))
colnames(sumNO3) <- c("Date","Treatment", "Applications", "meanNO3","seNO3")
pd <- position_dodge(0.3)
ggplot() +
#  geom_point(data = d, aes(x = Date, y = Nitrate_ppm, colour = Treatment),
#             size = 3, alpha = 0.2, position = pd) +
  geom_errorbar(data = sumNO3, aes(x = Date, ymin=meanNO3-seNO3, ymax=meanNO3+seNO3,
                                   colour = Applications),width=.4, position=pd, size =0.9) +
  geom_point(data = sumNO3, aes(Date, meanNO3, colour = Applications), 
             size = 4, position = pd) +
  geom_line(data = sumNO3, aes(Date, meanNO3, group = Treatment), colour = "black", 
            size = 1, position = pd) +
  labs(y = "[Nitrate] (ppm)", x = "Sample Date") +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(face="bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))



#### Bare Ground ####
# Note the datasheet has multiple sheets for different dates!
library(readxl)
library(dplyr)
rt <- read_xlsx("Chatfield_Row treatments.xlsx")
p <- read_xlsx("Veg_surveys_herbicide.xlsx", sheet = 2) %>%
  left_join(., rt, by = "Row") %>%
  filter(SPP == "G") %>%
  select(-Species, - `Seed Trt`) %>%
  group_by(Row, `Herb Trt`) %>%
  summarize(MeanCov = mean(Cover)) %>%
  mutate(Treatment = as.integer(`Herb Trt`)) %>%
  mutate(Treatment = as.character(Treatment)) %>%
  mutate(Treatment = as.factor(Treatment))

p <- read_xlsx("Veg_surveys_herbicide.xlsx") %>%
  left_join(., rt, by = "Row") %>%
  filter(SPP == "G") %>%
  select(-Species, - `Seed Trt`) %>%
  mutate(Treatment = as.integer(`Herb Trt`)) %>%
  mutate(Treatment = as.character(Treatment)) %>%
  mutate(Treatment = as.factor(Treatment))

leveneTest(p$MeanCov ~ p$Treatment)
m <- aov(p$MeanCov ~ p$Treatment)
shapiro.test(m$residuals)
summary(m) #
TukeyHSD(m)

ggplot(p, aes(Treatment, MeanCov)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5) +
  labs(x = "Herbicide Applications", 
       y = "% Bare Ground") +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))

p2 <- read_xlsx("Veg_surveys_herbicide.xlsx") %>%
  left_join(., rt, by = "Row") %>%
  filter(SPP != "G") %>%
  select(-Species, - `Seed Trt`) %>%
  group_by(Row, Quad, `Herb Trt`) %>%
  summarize(PlantCov = sum(Cover)) %>%
  ungroup() %>%
  group_by(Row, `Herb Trt`) %>%
  summarize(MeanCov = mean(PlantCov)) %>%
  mutate(Treatment = as.integer(`Herb Trt`)) %>%
  mutate(Treatment = as.character(Treatment)) %>%
  mutate(Treatment = as.factor(Treatment))

ggplot(p2, aes(Treatment, MeanCov)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5) +
  labs(x = "Herbicide Applications", 
       y = "% Plant Cover") +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        plot.margin = unit(c(0.1,0.1,0,0.1),"cm"))
