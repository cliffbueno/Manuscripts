# Analyzing Niwot Enzyme Data
setwd("~/Desktop/CU/1Classes/Arctic Alpine Biogeochem")
library(ggplot2)
library(plotrix)
library(pastecs)
library(WRS2)
library(WRS)
library(car)
source("~/Desktop/Functions/Summary.R")
data <- read.csv("NWT_Paired_Enzymes.csv")
mean(data$Dpinorg)
std.error(data$Dpinorg)
data <- subset(data, TIN != "NA")
mean(data$TIN)
std.error(data$TIN)

# 2002 Veg vs Unveg
two <- read.csv("2002.csv")
leveneTest(NL.PHOS ~ Treatment, data = two)
shapiro.test(two$NL.PHOS)
hist(two$NL.PHOS)
t.test(NL.PHOS ~ Treatment, data = two)
wilcox.test(NL.PHOS ~ Treatment, data = two)

leveneTest(log(NL.PHOS) ~ Treatment, data = two)
shapiro.test(log(two$NL.PHOS))
hist(log(two$NL.PHOS))
t.test(log(NL.PHOS) ~ Treatment, data = two)

leveneTest(BG.NL ~ Treatment, data = two)
shapiro.test(two$BG.NL)
hist(two$BG.NL)
t.test(BG.NL ~ Treatment, data = two)
wilcox.test(NL.PHOS ~ Treatment, data = two)

leveneTest(log(BG.NL) ~ Treatment, data = two)
shapiro.test(log(two$BG.NL))
hist(log(two$BG.NL))
t.test(log(BG.NL) ~ Treatment, data = two)

# Under 100 Stems. n = 65
data <- read.csv("NWT_Paired_Enzymes.csv")
sparse <- subset(data, Stems07 < 100)
sparse <- subset(sparse, Stems15 < 100) 
stat.desc(sparse, basic = FALSE, norm = TRUE)
sparse.nop <- subset(sparse, Sample.ID != 62)
sparse.non <- subset(sparse, Sample.ID != 78)
stat.desc(sparse.nop, basic = FALSE, norm = TRUE)
stat.desc(sparse.non, basic = FALSE, norm = TRUE)

m1 <- t.test(sparse$BG.NAG07, sparse$BG.NAG15, paired = TRUE, conf.level = 0.95) # Marginal p = 0.08644
m1
t1 <- m1$statistic[[1]]
df1 <- m1$parameter[[1]]
r1 <- sqrt(t1^2/(t1^2+df1))
round(r1, 3) # 0.22
wilcox.test(sparse$BG.NAG07, sparse$BG.NAG15, paired = TRUE) # Significant p = 0.03
mean(sparse$BG.NAG07, na.rm = TRUE)
mean(sparse$BG.NAG15, na.rm = TRUE) # Marginally higher
dif1 <- sparse$BG.NAG15 - sparse$BG.NAG07
shapiro.test(dif1) # Not normal
yuend(sparse$BG.NAG07, sparse$BG.NAG15) # p = 0.09
ydbt(sparse$BG.NAG07, sparse$BG.NAG15, nboot = 2000) # p = 0.08
bootdpci(sparse$BG.NAG07, sparse$BG.NAG15, est = tmean, nboot = 2000) # p = 0.03

m2 <- t.test(sparse$BG.PHOS07, sparse$BG.PHOS15, paired = TRUE, conf.level = 0.95) # NSD p = 0.31
m2
t2 <- m2$statistic[[1]]
df2 <- m2$parameter[[1]]
r2 <- sqrt(t2^2/(t2^2+df2))
round(r2, 3) # 0.137
wilcox.test(sparse$BG.PHOS07, sparse$BG.PHOS15, paired = TRUE) # NSD p = 0.64
mean(sparse$BG.PHOS07, na.rm = TRUE)
mean(sparse$BG.PHOS15, na.rm = TRUE)
dif2 <- sparse$BG.PHOS15 - sparse$BG.PHOS07
shapiro.test(dif2) # Not normal
yuend(sparse$BG.PHOS07, sparse$BG.PHOS15) # p = 0.62
ydbt(sparse$BG.PHOS07, sparse$BG.PHOS15, nboot = 2000) # p = 0.62
bootdpci(sparse$BG.PHOS07, sparse$BG.PHOS15, est = tmean, nboot = 2000) # p = 0.81

m3 <- t.test(sparse$NAG.PHOS07, sparse$NAG.PHOS15, paired = TRUE, conf.level = 0.95) # Marginal p = 0.05079
m3
t3 <- m3$statistic[[1]]
df3 <- m3$parameter[[1]]
r3 <- sqrt(t3^2/(t3^2+df3))
round(r3, 3) # 0.256
wilcox.test(sparse$NAG.PHOS07, sparse$NAG.PHOS15, paired = TRUE) # NSD p = 0.15
mean(sparse$NAG.PHOS07, na.rm = TRUE) # Marginally higher
mean(sparse$NAG.PHOS15, na.rm = TRUE)
dif3 <- sparse$BG.NAG15 - sparse$BG.NAG07
shapiro.test(dif3) # Not normal
yuend(sparse$NAG.PHOS07, sparse$NAG.PHOS15) # p = 0.04
ydbt(sparse$NAG.PHOS07, sparse$NAG.PHOS15, nboot = 2000) # p = 0.04
bootdpci(sparse$NAG.PHOS07, sparse$NAG.PHOS15, est = tmean, nboot = 2000) # p = 0.14

g <- read.csv("NWT_Paired_Enzymes_Sparse_Long.csv")
g <- subset(g, Ratio != "NA")
g$Year <- as.factor(g$Year)
sbgnag <- subset(g, Type == "BG:NAG")
leveneTest(Ratio ~ Year, data = sbgnag) # Good
sbgphos <- subset(g, Type == "BG:PHOS")
leveneTest(Ratio ~ Year, data = sbgphos) # 0.04, fine
snagphos <- subset(g, Type == "NAG:PHOS")
leveneTest(Ratio ~ Year, data = snagphos) # Good
sum <- summarySE(g, measurevar ="Ratio",groupvars=c("Year","Type"))
ggplot(sum, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_errorbar(data=sum,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  ylim(0.09,0.81) +
  theme_bw() +
  facet_wrap(~ Type, nrow = 3, scales = "free_y") +
  theme(strip.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16))

# New Graph with References
g <- read.csv("NWT_Paired_Enzymes_Sparse_Long_Ref.csv")
g <- subset(g, Ratio != "NA")
g$Year <- as.factor(g$Year)
l <- c('BG:NAG' = " ",'BG:PHOS' = " ", 'NAG:PHOS' = " ")
sum <- summarySE(g, measurevar ="Ratio",groupvars=c("Year","Type"))
ggplot(sum, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_errorbar(data=sum,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Type, nrow = 3, scales = "free_y", labeller = as_labeller(l)) +
  theme(strip.text = element_blank(),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        strip.background = element_blank())

############################## Over 100 Stems. n = 9 ########################################
data <- read.csv("NWT_Paired_Enzymes.csv")
dense <- subset(data, Stems07 > 100)
dense <- subset(dense, Stems15 > 100)
stat.desc(dense, basic = FALSE, norm = TRUE)
t.test(dense$BG.NAG07, dense$BG.NAG15, paired = TRUE, conf.level = 0.95) # Sig. p = 0.04472
mean(dense$BG.NAG07, na.rm = TRUE) # Significantly higher
mean(dense$BG.NAG15, na.rm = TRUE)
bootdpci(dense$BG.NAG07, dense$BG.NAG15, est = tmean, nboot = 2000) # p = 0.02

t.test(dense$BG.PHOS07, dense$BG.PHOS15, paired = TRUE, conf.level = 0.95) # NSD p = 0.6
mean(dense$BG.PHOS07, na.rm = TRUE)
mean(dense$BG.PHOS15, na.rm = TRUE)
bootdpci(dense$BG.PHOS07, dense$BG.PHOS15, est = tmean, nboot = 2000) # p = 0.53

t.test(dense$NAG.PHOS07, dense$NAG.PHOS15, paired = TRUE, conf.level = 0.95) # Marginal p = 0.06
mean(dense$NAG.PHOS07, na.rm = TRUE)
mean(dense$NAG.PHOS15, na.rm = TRUE) # Marginally higher
bootdpci(dense$NAG.PHOS07, dense$NAG.PHOS15, est = tmean, nboot = 2000) # p = 0

g2 <- read.csv("NWT_Paired_Enzymes_Dense_Long.csv")
g2 <- subset(g2, Ratio != "NA")
g2$Year <- as.factor(g2$Year)
dbgnag <- subset(g2, Type == "BG:NAG")
leveneTest(Ratio ~ Year, data = dbgnag) # Good
dbgphos <- subset(g2, Type == "BG:PHOS")
leveneTest(Ratio ~ Year, data = dbgphos) # Good
dnagphos <- subset(g2, Type == "NAG:PHOS")
leveneTest(Ratio ~ Year, data = dnagphos) # Good
sum <- summarySE(g2, measurevar ="Ratio",groupvars=c("Year","Type"))
ggplot(sum, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_errorbar(data=sum,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Type, nrow = 3, scales = "free_y") +
  theme(strip.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16))



######################### Plots that were under 100 now over 100. n = 21 #######
data <- read.csv("NWT_Paired_Enzymes.csv")
change <- subset(data, Stems07 < 100)
change <- subset(change, Stems15 > 100)
stemdiff <- change$Stems15 - change$Stems07
mean(stemdiff)
std.error(stemdiff)
t.test(change$BG.NAG07, change$BG.NAG15, paired = TRUE, conf.level = 0.95) # NSD 0.81
mean(change$BG.NAG07, na.rm = TRUE)
mean(change$BG.NAG15, na.rm = TRUE)
bootdpci(change$BG.NAG07, change$BG.NAG15, est = tmean, nboot = 2000) # p = 0.95

t.test(change$BG.PHOS07, change$BG.PHOS15, paired = TRUE, conf.level = 0.95) # NSD p = 0.23
mean(change$BG.PHOS07, na.rm = TRUE)
mean(change$BG.PHOS15, na.rm = TRUE)
bootdpci(change$BG.PHOS07, change$BG.PHOS15, est = tmean, nboot = 2000) # p = 0.24

t.test(change$NAG.PHOS07, change$NAG.PHOS15, paired = TRUE, conf.level = 0.95) # NSD p = 0.2092
mean(change$NAG.PHOS07, na.rm = TRUE)
mean(change$NAG.PHOS15, na.rm = TRUE)
bootdpci(change$NAG.PHOS07, change$NAG.PHOS15, est = tmean, nboot = 2000) # p = 0.32

mean(change$Stems15 - change$Stems07)

g3 <- read.csv("NWT_Paired_Enzymes_Change_Long.csv")
g3 <- subset(g3, Ratio != "NA")
g3$Year <- as.factor(g3$Year)
cbgnag <- subset(g3, Type == "BG:NAG")
leveneTest(Ratio ~ Year, data = cbgnag) # Good
cbgphos <- subset(g3, Type == "BG:PHOS")
leveneTest(Ratio ~ Year, data = cbgphos) # 0.04, fine
cnagphos <- subset(g3, Type == "NAG:PHOS")
leveneTest(Ratio ~ Year, data = cnagphos) # Good
sum <- summarySE(g3, measurevar ="Ratio",groupvars=c("Year","Type"))
ggplot(sum, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_errorbar(data=sum,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Type, nrow = 3, scales = "free_y") +
  theme(strip.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16))

############################ Continuous ####################################
d <- read.csv("NWT_Paired_Enzymes.csv")

m1 <- lm(BG.NAG15 ~ Stems15, data = d)
summary(m1) # Sig negative, p = 0.04, R2 = 0.05
plot(BG.NAG15 ~ Stems15, data = d)
abline(m1)

m2 <- lm(BG.PHOS15 ~ Stems15, data = d)
summary(m2) # NR p = 0.9
plot(BG.PHOS15 ~ Stems15, data = d)
abline(m2)

m3 <- lm(NAG.PHOS15 ~ Stems15, data = d)
summary(m3) # Sig positive, p = 0.007, R2 = 0.08
plot(NAG.PHOS15 ~ Stems15, data = d)
abline(m3)



#################### New Multipanel Graphs by Ratio and Veg ###################
v <- read.csv("BG.NAG.csv")
v <- subset(v, Ratio != "NA")
v$Veg <- factor(v$Veg, levels = c("Sparse","Change","Dense"))
v$Year <- as.factor(v$Year)
sum <- summarySE(v, measurevar ="Ratio",groupvars=c("Veg","Year"))
ggplot(sum, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 1.93, linetype = "dotted") +
  geom_hline(yintercept = 1.827, linetype = "longdash") +
  geom_hline(yintercept = 7.5, linetype = "dotdash") +
  geom_errorbar(data=sum,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Veg, ncol = 3) +
  theme(strip.text = element_blank(),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        strip.background = element_blank())

v1 <- read.csv("BG.PHOS.csv")
v1 <- subset(v1, Ratio != "NA")
v1$Veg <- factor(v1$Veg, levels = c("Sparse","Change","Dense"))
v1$Year <- as.factor(v1$Year)
sum1 <- summarySE(v1, measurevar ="Ratio",groupvars=c("Veg","Year"))
ggplot(sum1, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.62, linetype = "dotted") +
  geom_hline(yintercept = 0.214, linetype = "longdash") +
  geom_hline(yintercept = 0.845, linetype = "dotdash") +
  geom_errorbar(data=sum1,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Veg, ncol = 3) +
  theme(strip.text = element_blank(),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        strip.background = element_blank())


v2 <- read.csv("NAG.PHOS.csv")
v2 <- subset(v2, Ratio != "NA")
v2$Veg <- factor(v2$Veg, levels = c("Sparse","Change","Dense"))
v2$Year <- as.factor(v2$Year)
sum2 <- summarySE(v2, measurevar ="Ratio",groupvars=c("Veg","Year"))
ggplot(sum2, aes(Year, Ratio)) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.32, linetype = "dotted") +
  geom_hline(yintercept = 0.126, linetype = "longdash") +
  geom_hline(yintercept = 0.175, linetype = "dotdash") +
  geom_errorbar(data=sum2,aes(ymin=Ratio-se,ymax=Ratio+se),width=.3,size=0.9) +
  labs(x = "Year", y = "Enzyme Ratio") +
  theme_bw() +
  facet_wrap(~ Veg, ncol = 3) +
  theme(strip.text = element_blank(),
        axis.title = element_text(face="bold", size = 18), 
        axis.text = element_text(size = 16),
        panel.grid = element_blank(),
        strip.background = element_blank())
