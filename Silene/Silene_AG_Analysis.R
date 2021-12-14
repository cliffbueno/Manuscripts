# Analysis of Silene acaulis growth and fungal colonization
# Cliff Bueno de Mesquita and Conor Meade, July 2017
# Irish and Colorado Silene grown in alpine room in Irish and Colorado soils March - July

# Load packages
library(plyr)
library(tidyverse)
library(PMCMR)
library(car)
library(reshape2)

# Load data
setwd("~/Desktop/CU/2Research/Silene")
d <- read.csv("Silene_Data.csv")
scope <- subset(d, Microscopy == "Yes")
seq <- subset(d, SequenceData == "Yes")
table(d$Treatment)
mean(d$Days)
(sd(d$Days))/(sqrt(nrow(d)))

############################ Above Ground Plant Growth ################################
# Biomass, Area, Height, Longest Leaf, Leaf Number
# Make normal, then run anova tests and graph

# Explore data
hist(d$Days)

plot(d$Treatment, d$Biomass)
plot(d$Treatment, d$Area)
plot(d$Treatment, d$Height)
plot(d$Treatment, d$Longest)
plot(d$Treatment, d$Leafnum)

plot(d$Treatment, d$Growth)
plot(d$Treatment, d$A.Day)
plot(d$Treatment, d$H.Day)
plot(d$Treatment, d$LL.Day)
plot(d$Treatment, d$LN.Day)

# Check normality (none normal)
shapiro.test(d$Growth)
hist(d$Growth)
shapiro.test(d$A.Day)
hist(d$A.Day)
shapiro.test(d$H.Day)
hist(d$H.Day)
shapiro.test(d$LL.Day)
hist(d$LL.Day)
shapiro.test(d$LN.Day)
hist(d$LN.Day)

# Transform to at least approximate normality to use ANOVA test
shapiro.test(log(d$Growth))
hist(log(d$Growth)) # Fine
d$GrowthLog <- log(d$Growth)

shapiro.test(log(d$A.Day))
hist(log(d$A.Day)) # Great
d$A.Day.Log <- log(d$A.Day)

shapiro.test(log(d$H.Day))
hist(log(d$H.Day)) # Good
d$H.Day.Log <- log(d$H.Day)

shapiro.test(log(d$LL.Day))
hist(log(d$LL.Day)) # Good
d$LL.Day.Log <- log(d$LL.Day)

shapiro.test(log(d$LN.Day))
hist(log(d$LN.Day)) # Fine
d$LN.Day.Log <- log(d$LN.Day)


# BIOMASS per day
m <- aov(d$GrowthLog ~ d$Plant*d$Soil)
summary(m) # Significant soil effect and interaction, marginal plant effect
Anova(m)
TukeyHSD(m)
plot(m)

growth1 <- ddply(d, .(Plant, Soil), summarise, GrowthLog = mean(GrowthLog))
summ <- ddply(d, .(Plant, Soil), summarise, GrowthLog = mean(GrowthLog))
ggplot(d, aes(as.factor(Soil), GrowthLog, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("ln(Biomass(g)/Day)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.86),
        legend.background = element_rect(colour = "black"))

## AREA per day
m <- aov(d$A.Day.Log ~ d$Plant*d$Soil)
summary(m) # Significant plant effect and interaction
TukeyHSD(m)
plot(m)

growth2 <- ddply(d, .(Plant, Soil), summarise, GrowthLog = mean(A.Day.Log))
summ <- ddply(d, .(Plant, Soil), summarise, A.Day.Log = mean(A.Day.Log))
ggplot(d, aes(as.factor(Soil), A.Day.Log, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("ln(Area(mm^2)/Day)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.86),
        legend.background = element_rect(colour = "black"))

## HEIGHT per day
m <- aov(d$H.Day.Log ~ d$Plant*d$Soil)
summary(m) # Significant plant effect and interaction
TukeyHSD(m)
plot(m)

growth3 <- ddply(d, .(Plant, Soil), summarise, GrowthLog = mean(H.Day.Log))
summ <- ddply(d, .(Plant, Soil), summarise, H.Day.Log = mean(H.Day.Log))
ggplot(d, aes(as.factor(Soil), H.Day.Log, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("ln(Height(mm/Day)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.86),
        legend.background = element_rect(colour = "black"))

## LONGEST LEAF LENGTH per day
m <- aov(d$LL.Day.Log ~ d$Plant*d$Soil)
summary(m) # Significant plant effect and interaction
TukeyHSD(m)
plot(m)

summ <- ddply(d, .(Plant, Soil), summarise, LL.Day.Log = mean(LL.Day.Log))
ggplot(d, aes(as.factor(Soil), LL.Day.Log, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("ln(Longest Leaf Growth Rate (mm/Day))") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.86),
        legend.background = element_rect(colour = "black"))

## LEAF NUMBER per day
m <- aov(d$LN.Day.Log ~ d$Plant*d$Soil)
summary(m) # All significant
TukeyHSD(m)
plot(m)

growth4 <- ddply(d, .(Plant, Soil), summarise, GrowthLog = mean(LN.Day.Log))
summ <- ddply(d, .(Plant, Soil), summarise, LN.Day.Log = mean(LN.Day.Log))
ggplot(d, aes(as.factor(Soil), LN.Day.Log, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("ln(New Leaves Per Day)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.86),
        legend.background = element_rect(colour = "black"))

# Multipanel plant response graph
# If means needed
growth1$dataset <- "a) Biomass"
growth2$dataset <- "b) Area"
growth3$dataset <- "c) Height"
growth4$dataset <- "d) Leaf Number"
growth_mp <- rbind(growth1, growth2, growth3, growth4)

# Melt by the four response variables, then facet_wrap
d_melt <- melt(d, 
               measure.vars = c("GrowthLog", "A.Day.Log", "H.Day.Log", "LN.Day.Log"),
               id.vars = c("Soil", "Plant")) %>%
  mutate(variable = revalue(variable,
          c("GrowthLog" = "a) Log Biomass (g/day)",
          "A.Day.Log" = "b) Log Area (mm/day)",
          "H.Day.Log" = "c) Log Height (mm/day)",
          "LN.Day.Log" = "d) Log Leaf Number (number/day)")))
pdf("GrowthMultipanel.pdf", width = 7, height = 4.5)
ggplot(d_melt, aes(Plant, value, colour = Soil)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.25, position = position_jitterdodge()) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  labs(x = "Plant", y = NULL) +
  facet_wrap( ~ variable, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 11),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))
dev.off()



################################# Bray - Curtis Multipanel #################################
fbc <- read.csv("ITS_Dissimilarity.csv")
fbc$dataset <- "a) Fungi"
bbc <- read.csv("16S_Dissimilarity.csv")
bbc$dataset <- "b) Bacteria"
bc <- rbind(fbc, bbc)
bc$Comparison = factor(bc$Comparison, levels=c("COsoilCOfield", "COsoilCOgh", "COsoilIRgh",
                                               "IRsoilIRfield", "IRsoilIRgh","IRsoilCOgh"))
toremove <- which(bc$Comparison=="COsoilCOfield"|bc$Comparison=="IRsoilIRfield")
bc <- droplevels(bc[-toremove,])
pdf("BCMultipanel.pdf", width = 6, height = 3)
ggplot(bc, aes(Plant, Dissimilarity, colour = Soil)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, alpha = 0.25, position = position_jitterdodge()) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  labs(x = "Plant", y = "Bray-Curtis Dissimilarity") +
  facet_wrap( ~ dataset, ncol = 2) +
  ylim(0.5, 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 11),
        legend.position = c(0.1, 0.225),
        legend.background = element_rect(colour = "black"))
dev.off()



################################### Root Colonization ########################################
plot(scope$Treatment, scope$DSE)
shapiro.test(scope$DSE)
shapiro.test(log(scope$DSE))
hist(scope$DSE)
hist(log(scope$DSE))
d$GrowthLog <- log(d$Growth)
m <- aov(scope$DSE ~ scope$Plant*scope$Soil)
summary(m) # Irish plant marginally higher
TukeyHSD(m)
plot(scope$Soil, scope$DSE)
plot(scope$Plant, scope$DSE)

summ <- ddply(scope, .(Plant, Soil), summarise, DSE = mean(DSE))
ggplot(scope, aes(as.factor(Soil), DSE, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("% Colonization") +
  ylim(0,20) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        legend.position = c(0.91,0.88),
        legend.background = element_rect(colour = "black"))

plot(scope$Treatment, scope$AMF)
shapiro.test(scope$AMF)
shapiro.test(log(scope$AMF + 1))
hist(scope$AMF)
hist(log(scope$AMF + 1))
scope$logAMF <- log(scope$AMF + 1)
m <- aov(scope$logAMF ~ scope$Plant*scope$Soil)
summary(m) # NSD
TukeyHSD(m)
plot(scope$Soil, scope$AMF)
plot(scope$Plant, scope$AMF)
plot(scope$DSE, scope$Growth)
plot(scope$AMF, scope$Growth)

summ <- ddply(scope, .(Plant, Soil), summarise, AMF = mean(AMF))
ggplot(scope, aes(as.factor(Soil), AMF, colour = Plant)) +
  geom_boxplot() +
  geom_point(data = summ, aes(group=Plant), colour = "black",
             position = position_dodge(width=0.75)) +
  geom_line(data = summ, aes(group=Plant),
            position = position_dodge(width=0.75)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Soil") +
  ylab("% Colonization") +
  ylim(0,20) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# New Graph
ggplot(scope, aes(as.factor(Plant), AMF, colour = Soil)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Plant") +
  ylab("% Colonization") +
  theme_bw() +
  ylim(0,20) +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = c(0.9,0.88),
        legend.background = element_rect(colour = "black"))

ggplot(scope, aes(as.factor(Plant), DSE, colour = Soil)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_colour_manual(values = c("blue2","darkorange2")) +
  xlab("Plant") +
  ylab("% Colonization") +
  theme_bw() +
  ylim(0,20) +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        legend.position = c(0.91,0.88),
        legend.background = element_rect(colour = "black"))



############################## Sequence Relative Abundances #####################################
seq <- na.omit(seq)
plot(seq$Mut.Abund, seq$GrowthLog)
plot(seq$Path.Abund, seq$GrowthLog)
plot(seq$DSE, seq$GrowthLog)
plot(seq$AMF, seq$GrowthLog)
