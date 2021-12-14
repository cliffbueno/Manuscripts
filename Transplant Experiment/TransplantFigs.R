# Cliff transplant paper biomass and survival figures
library(tidyverse)
library(scales)
library(devtools)
source_url("https://github.com/cliffbueno/Functions/blob/master/Summary.R?raw=TRUE")

#### Combined Multipanel Biomass Figure ####
setwd("~/Desktop/CU/2Research/Transplant/Field/Year2")
alldes <- read.csv("Des_2016to2018_Pots.csv")
alldes$Year <- as.factor(alldes$Year)
alldes$Block <- as.factor(alldes$Block)
desgh <- subset(alldes, Year == "2016")
des17 <- subset(alldes, Year == "2017")
des18 <- subset(alldes, Year == "2018")
field <- subset(alldes, Year != "2016")
desgh$Ind <- as.factor(desgh$Ind)
des17$Ind <- as.factor(des17$Ind)
des18$Ind <- as.factor(des18$Ind)

# Combined Multipanel Biomass Figure
U1yr1 <- subset(des17, Inoculum == "BA")
U2yr1 <- subset(des17, Inoculum == "BB")
U3yr1 <- subset(des17, Inoculum == "BC")
U4yr1 <- subset(des17, Inoculum == "BD")
V1yr1 <- subset(des17, Inoculum == "TA")
V2yr1 <- subset(des17, Inoculum == "TB")
V3yr1 <- subset(des17, Inoculum == "TC")
V4yr1 <- subset(des17, Inoculum == "TD")
U1yr2 <- subset(des18, Inoculum == "BA")
U2yr2 <- subset(des18, Inoculum == "BB")
U3yr2 <- subset(des18, Inoculum == "BC")
U4yr2 <- subset(des18, Inoculum == "BD")
V1yr2 <- subset(des18, Inoculum == "TA")
V2yr2 <- subset(des18, Inoculum == "TB")
V3yr2 <- subset(des18, Inoculum == "TC")
V4yr2 <- subset(des18, Inoculum == "TD")
show_col(hue_pal()(9))

# New Biomass Facet Wrap, GSL on X, show interactions
facet_names <- c(`2017` = "a) Field Year 1",`2018` = "b) Field Year 2")
pd <- position_dodge(0.3)
ggplot(field,aes(GSL, BiomassPerIndDiff)) + 
  geom_point(size = 4, aes(colour = Inoculum, shape = Snow), position = pd) + 
  geom_smooth(data = U1yr1, method = lm, se = FALSE, col = "#F8766D", size = 0.5) +
  geom_smooth(data = U1yr2, method = lm, se = FALSE, col = "#F8766D", size = 0.5) +
  geom_smooth(data = U2yr1, method = lm, se = FALSE, col = "#CD9600", size = 0.5) +
  geom_smooth(data = U2yr2, method = lm, se = FALSE, col = "#CD9600", size = 0.5) +
  geom_smooth(data = U3yr1, method = lm, se = FALSE, col = "#7CAE00", size = 0.5) +
  geom_smooth(data = U3yr2, method = lm, se = FALSE, col = "#7CAE00", size = 0.5) +
  geom_smooth(data = U4yr1, method = lm, se = FALSE, col = "#00BE67", size = 0.5) +
  geom_smooth(data = U4yr2, method = lm, se = FALSE, col = "#00BE67", size = 0.5) +
  geom_smooth(data = V1yr1, method = lm, se = FALSE, col = "#00BFC4", size = 0.5) +
  geom_smooth(data = V1yr2, method = lm, se = FALSE, col = "#00BFC4", size = 0.5) +
  geom_smooth(data = V2yr1, method = lm, se = FALSE, col = "#00A9FF", size = 0.5) +
  geom_smooth(data = V2yr2, method = lm, se = FALSE, col = "#00A9FF", size = 0.5) +
  geom_smooth(data = V3yr1, method = lm, se = FALSE, col = "#C77CFF", size = 0.5) +
  geom_smooth(data = V3yr2, method = lm, se = FALSE, col = "#C77CFF", size = 0.5) +
  geom_smooth(data = V4yr1, method = lm, se = FALSE, col = "#FF61CC", size = 0.5) +
  geom_smooth(data = V4yr2, method = lm, se = FALSE, col = "#FF61CC", size = 0.5) +
  scale_colour_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF","#FF61CC"), labels = c("BA"="U1","BB"="U2","BC"="U3","BD"="U4","TA"="V1","TB"="V2","TC"="V3","TD"="V4")) +
  scale_shape_manual(name = "Plot", values = c(16,17)) +
  labs(y = "Growth (g/ind./yr)", x = "Growing Season (Days)") +
  theme_bw() +
  facet_wrap(~Year, nrow = 2, scales = "free_y") +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(face="bold", size = 18),
        legend.text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

pd <- position_dodge(0.5)
sum<-summarySE(field,measurevar="BiomassPerIndDiff",groupvars=c("Year","Treatment","Inoculum"))
sum$I <- c("I","I","I","I","I","I","I","I","I","I","I","I","I","I","I","I")
ggplot(sum, aes(I, BiomassPerIndDiff, colour = Inoculum)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(data=sum,aes(ymin=BiomassPerIndDiff-se,ymax=BiomassPerIndDiff+se),width=.05, position = pd) +
  labs(x = "I", y = "Biomass (g/individual/yr)") +
  theme_bw() +
  facet_wrap(~Year, scales = "free_y", ncol = 1, nrow = 2) +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18),
        strip.text = element_text(size = 16))



#### Combined Multipanel Survival Figure ####
setwd("~/Desktop/CU/2Research/Transplant/Field/")
ppc <- read.csv("PredProbCombinedglm.csv")
ppc$Year <- as.factor(ppc$Year)
pp17 <- subset(ppc, Year == "2017")
pp18 <- subset(ppc, Year == "2018")
U1yr1 <- subset(pp17, Inoculum == "U1")
U2yr1 <- subset(pp17, Inoculum == "U2")
U3yr1 <- subset(pp17, Inoculum == "U3")
U4yr1 <- subset(pp17, Inoculum == "U4")
V1yr1 <- subset(pp17, Inoculum == "V1")
V2yr1 <- subset(pp17, Inoculum == "V2")
V3yr1 <- subset(pp17, Inoculum == "V3")
V4yr1 <- subset(pp17, Inoculum == "V4")
U1yr2 <- subset(pp18, Inoculum == "U1")
U2yr2 <- subset(pp18, Inoculum == "U2")
U3yr2 <- subset(pp18, Inoculum == "U3")
U4yr2 <- subset(pp18, Inoculum == "U4")
V1yr2 <- subset(pp18, Inoculum == "V1")
V2yr2 <- subset(pp18, Inoculum == "V2")
V3yr2 <- subset(pp18, Inoculum == "V3")
V4yr2 <- subset(pp18, Inoculum == "V4")
# show_col(hue_pal()(8))

# New Survival Facet Wrap, GSL on X, show interactions
pd <- position_dodge(0.3)
ggplot(ppc,aes(GSL, PP)) + 
  geom_point(size = 4, aes(colour = Inoculum, shape = Snow), position = pd) + 
  geom_smooth(data = U1yr1, method = lm, se = FALSE, col = "#F8766D", size = 0.5) +
  geom_smooth(data = U1yr2, method = lm, se = FALSE, col = "#F8766D", size = 0.5) +
  geom_smooth(data = U2yr1, method = lm, se = FALSE, col = "#CD9600", size = 0.5) +
  geom_smooth(data = U2yr2, method = lm, se = FALSE, col = "#CD9600", size = 0.5) +
  geom_smooth(data = U3yr1, method = lm, se = FALSE, col = "#7CAE00", size = 0.5) +
  geom_smooth(data = U3yr2, method = lm, se = FALSE, col = "#7CAE00", size = 0.5) +
  geom_smooth(data = U4yr1, method = lm, se = FALSE, col = "#00BE67", size = 0.5) +
  geom_smooth(data = U4yr2, method = lm, se = FALSE, col = "#00BE67", size = 0.5) +
  geom_smooth(data = V1yr1, method = lm, se = FALSE, col = "#00BFC4", size = 0.5) +
  geom_smooth(data = V1yr2, method = lm, se = FALSE, col = "#00BFC4", size = 0.5) +
  geom_smooth(data = V2yr1, method = lm, se = FALSE, col = "#00A9FF", size = 0.5) +
  geom_smooth(data = V2yr2, method = lm, se = FALSE, col = "#00A9FF", size = 0.5) +
  geom_smooth(data = V3yr1, method = lm, se = FALSE, col = "#C77CFF", size = 0.5) +
  geom_smooth(data = V3yr2, method = lm, se = FALSE, col = "#C77CFF", size = 0.5) +
  geom_smooth(data = V4yr1, method = lm, se = FALSE, col = "#FF61CC", size = 0.5) +
  geom_smooth(data = V4yr2, method = lm, se = FALSE, col = "#FF61CC", size = 0.5) +
  scale_shape_manual(name = "Plot", values = c(16,17)) +
  ylim(0.1, 1) +
  labs(y = "Predicted Probability of Survival", x = "Growing Season (Days)") +
  theme_bw() +
  facet_wrap(~Year, nrow = 2, labeller = as_labeller(facet_names)) +
  theme(legend.position = "none",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(face="bold", size = 18),
        legend.text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18))

pd <- position_dodge(0.5)
sum<-summarySE(ppc,measurevar="PP",groupvars=c("Year","Treatment","Inoculum"))
sum$I <- c("I","I","I","I","I","I","I","I","I","I","I","I","I","I","I","I")
ggplot(sum, aes(I, PP,colour = Inoculum)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(data=sum,aes(ymin=PP-se,ymax=PP+se),width=.1,position=pd) +
  labs(x = "I", y = "Predicted Probability of Survival") +
  theme_bw() +
  ylim(0.1,1) +
  facet_wrap(~Year, scales = "free_y", ncol = 1, nrow = 2) +
  theme(legend.position = "right",
        legend.background = element_rect(colour = "black"),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(face="bold", size = 18), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(face="bold",size=18),
        strip.text = element_text(size = 16))