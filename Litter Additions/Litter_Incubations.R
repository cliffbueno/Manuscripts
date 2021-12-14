# Analyze litter lab incubation experiment
# 5 jars of des, oxy, sil, control.
# soil moistened to 90% WHC
# incubated at 10˚C
# measured with EGM 4 twice a week

# Get Started
library(Matrix)
library(lme4)
library(plyr)
library(ggplot2)
library(nlme)
library(car)
library(rcompanion)
library(multcomp)
library(scales)
library(pbkrtest)
library(Rmixmod)
source("~/Desktop/Functions/Summary.R")
setwd("~/Desktop/CU/2Research/Litter/Lab")
d <- read.csv("Litter_Incubations.csv")
cont <- subset(d, Treatment == "Control")
oxy <- subset(d, Treatment == "Oxyria")
des <- subset(d, Treatment == "Deschampsia")
sil <- subset(d, Treatment == "Silene")

# Test first vs second weekly sample
buildup <- subset(d, WeekPoint == "1" | WeekPoint == "2")
buildup$WeekPoint <- as.factor(buildup$WeekPoint)
summary(aov(buildup$WeekPointBuildup ~ buildup$Treatment+buildup$WeekPoint))
summary(aov(buildup$DailyBuildup ~ buildup$Treatment+buildup$WeekPoint))
t.test(buildup$DailyBuildup ~ buildup$WeekPoint)

m <- aov(CO2_Cum ~ Days*Treatment, data = d)
summary(m)

show_col(hue_pal()(4))

# Repeated Measures ANOVA
model <- gls(CO2_Cum ~ Treatment + Days + Treatment*Days, correlation = corAR1(form = ~ Days | JarID),data=d,method="REML")
Anova(model) # Significant Treatment, Days, Interaction
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}

multCompTukey <- glht(model, linfct = mcp(Treatment = "Tukey"))
summary(multCompTukey)

final <- d[201:220,]
m <- aov(CO2_Cum ~ Treatment, data = final)
summary(m)
TukeyHSD(m)

summary(lm(oxy$CO2_Cum ~ poly(oxy$Days,2,raw=TRUE)))
o <- function (x) -69*x^2 + 4369.75*x + 1915.54

summary(lm(des$CO2_Cum ~ poly(des$Days,2,raw=TRUE)))
de <- function (x) -20.23*x^2 + 1862.59*x + 2700.67

summary(lm(sil$CO2_Cum ~ poly(sil$Days,2,raw=TRUE)))
s <- function (x) -13.95*x^2 + 1769.95*x - 527.42

summary(lm(cont$CO2_Cum ~ poly(cont$Days,2,raw=TRUE)))
c <- function (x) 0.02198*x^2 + 76.1498*x + 478.5323

d$Days <- as.factor(d$Days)
sum <- summarySE(d, measurevar = "CO2_Cum", groupvars = "Days")
pd <- position_dodge(1.5)
ggplot(sum, aes(Days, CO2_Cum, shape = Treatment, group = Treatment)) +
  stat_function(fun=o, colour="black") +
  stat_function(fun=de, colour="black") +
  stat_function(fun=s, colour="black") +
  stat_function(fun=c, colour="black") +
  geom_errorbar(aes(ymin=CO2_Cum-se, ymax=CO2_Cum+se), width=.1, position=pd) +
  geom_point(position = pd, size = 3) +
  labs(x = "Day Incubated", y = expression(CO[2]*" (ppm)")) +
  scale_shape_manual(values = 4:7) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold",size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(face="bold", size=16))


############################## Analyze Daily Rate from each week ################################
d <- subset(d, Rate != "NA")
d <- subset(d, Rate != 0)
cont <- subset(d, Treatment == "Control")
oxy <- subset(d, Treatment == "Oxyria")
des <- subset(d, Treatment == "Deschampsia")
sil <- subset(d, Treatment == "Silene")
show_col(hue_pal()(4))

m <- aov(Rate ~ Days*Treatment, data = d)
summary(m) # All significant
shapiro.test(m$residuals)

# Repeated Measures ANOVA
model <- gls(Rate ~ Treatment + Days + Treatment*Days, correlation = corAR1(form = ~ Days | JarID), data=d , method="REML")
Anova(model) # Significant Treatment, Days, Interaction
Anova(model, test.statistic = "F") # Doesn't work
anovaTAB(model,d)
x = residuals(model)
plotNormalHistogram(x)
plot(fitted(model),residuals(model))

model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)
}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)
}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)
}

multCompTukey <- glht(model, linfct = mcp(Treatment = "Tukey"))
summary(multCompTukey)

final <- subset(d, Days == 35)
m <- aov(Rate ~ Treatment, data = final)
summary(m)
TukeyHSD(m)

# For graph include initial zero
d <- read.csv("Litter_Incubations.csv")
d <- subset(d, Rate != "NA")
cont <- subset(d, Treatment == "Control")
oxy <- subset(d, Treatment == "Oxyria")
des <- subset(d, Treatment == "Deschampsia")
sil <- subset(d, Treatment == "Silene")
show_col(hue_pal()(4))

summary(lm(oxy$Rate ~ poly(oxy$Days,4,raw=TRUE)))
o <- function (x) 1151.57203*x - 106.32670*x^2 + 3.36275*x^3 - 0.03556*x^4 - 6.58231

summary(lm(des$Rate ~ poly(des$Days,4,raw=TRUE)))
de <- function (x) 558.78891*x - 53.24676*x^2 + 1.79813*x^3 - 0.02040*x^4 + 13.42052

summary(lm(sil$Rate ~ poly(sil$Days,2,raw=TRUE)))
s <- function (x) 135.037*x - 3.444*x^2 + 201.993

summary(lm(cont$Rate ~ poly(cont$Days,2,raw=TRUE)))
c <- function (x) 6.55061*x - 0.12749*x^2 + 17.20408

sum <- summarySE(d, measurevar = "Rate", groupvars = c("Days", "Treatment"))
pd <- position_dodge(1)
ggplot(sum, aes(Days, Rate, colour = Treatment, group = Treatment)) +
  stat_function(fun=c, colour="#F8766D") +
  stat_function(fun=de, colour="#7CAE00") +
  stat_function(fun=o, colour="#00BFC4") +
  stat_function(fun=s, colour="#C77CFF") +
  geom_errorbar(aes(ymin=Rate-se, ymax=Rate+se), width=.1, position=pd) +
  geom_point(position = pd, size = 3) +
  labs(x = "Days Incubated", y = expression(bold(CO[2]*" (ppm/day)"))) +
  theme_bw() +
  theme(axis.title = element_text(face="bold",size = 18), 
        axis.text = element_text(size = 16))
