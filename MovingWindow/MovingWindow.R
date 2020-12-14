# NWT PCA Moving Window Analysis and Other PCA Analyses
# by Cliff Bueno de Mesquita, September/October 2019
# Goal: Test how the axes and loadings change over the course of 10 year windows
# There will be 28 windows analyzed (1982-1991 through 2009-2018)
# Over these 28 points, look at how axes scores of the variables change
# 19 years will have axes scores for 11 windows, look at how axes scores of those years change
# Also look at how the PCA performs (ie how much variation the first two axes explain)
# Also try removing variables, particulary Ice Off (other sites don't have), maybe others
# Try increasing window size to 15 yrs
# Get D1 winter precip data. Look for correlations over windows with temp and GSL variables
# Get D1 summer temp data and check correlations with that
# Reviewer response - update figures 2, 3, 5, 6. Plot MW correlations for all variables with temperature over time. Plot PC1 scores for all variables over time.



##################################### Data Setup ###############################
library(plyr)
library(tidyverse)
library(vegan)
library(reshape2)
library(gridExtra)
library(scales)
library(breakpoint)
library(FSA)
library(hornpa)
library(car)
library(segmented)
theme_set(theme_bw())
setwd("~/Documents/GitHub/long-term-trends/extended_summer/sensitivity_analyses")
NWTclimate <- read.csv("../NWT_sumallPCclimate_19822018.csv")
sumallyrsOutput <- NWTclimate[c("eco_year", "sumallPC1", "sumallPC2")]
sumallyrsVarout <- read.csv("../analysis/output_data/extsum_pca_input/crall21x/NWT_sumallPCvarout_19822018.csv")
sumallyrsVarout$variable
sumallyrsVarout$varshort <- c("Temp", "Precip", "MoistureDeficit", "PET", "GDD", 
                              "Days5C", "Days12C", "GSL", "IceOff")
sumallyrsVarout$varnum[sumallyrsVarout$variable == "sum_meanT"] <- 3
sumallyrsVarout$varnum[sumallyrsVarout$variable == "sum_precip"] <- 1
sumallyrsVarout$varnum[sumallyrsVarout$variable == "sum_moisturedeficit"] <- 6
sumallyrsVarout$varnum[sumallyrsVarout$variable == "sum_PET"] <- 4
sumallyrsVarout$varnum[sumallyrsVarout$variable == "sum_GDD"] <- 2
sumallyrsVarout$varnum[sumallyrsVarout$variable == "fivedayrunning5C"] <- 8
sumallyrsVarout$varnum[sumallyrsVarout$variable == "fivedayrunning12C"] <- 9
sumallyrsVarout$varnum[sumallyrsVarout$variable == "GSLthreedayneg3C"] <- 5
sumallyrsVarout$varnum[sumallyrsVarout$variable == "iceoff_GL4"] <- 7
sumallyrsVarout <- sumallyrsVarout %>%
  dplyr::rename(sumallPC1 = PC1,
                sumallPC2 = PC2)

# Subset climate data for summer variables
climateSummer <- NWTclimate %>%
  dplyr::rename(year = eco_year) %>%
  select(year, sum_meanT, sum_precip, sum_moisturedeficit, sum_PET, sum_GDD, 
         fivedayrunning5C, fivedayrunning12C, 
         GSLthreedayneg3C, iceoff_GL4) %>%
  na.omit() # no NAs allowed (there aren't any in ctw's updated 1982-2017 anyway)
row.names(climateSummer) <- climateSummer$year
climateSummer <- as.data.frame(climateSummer) # preserves rownames in PCA output
mean(climateSummer$sum_meanT)
se(climateSummer$sum_meanT)

# Full PCA
sumallPCA <- rda(na.exclude(climateSummer[,2:ncol(climateSummer)]), scale=T)
s <- summary(sumallPCA) 
s 
# PC1 (extended summer) explains 52.79%, PC2 (moisture gradient) explains 18.53%
par(mfrow = c(1,1))
biplot(sumallPCA, main = "1982-ongoing")

# Parallel Analysis for number of components to retain (Horn, 1965)
set.seed(921)
hornpa(k = 9, size = 37, reps = 1000, seed = 921)

# Calculate Correlations from Loadings 
loadings <- scores(sumallPCA, display = "species", scaling = 0)
evs <- eigenvals(sumallPCA)
evs # PC1 and PC2 are greater than mean values from parallel analysis!
climateSummerscaled <- as.data.frame(scale(climateSummer[,2:ncol(climateSummer)]))
# Axis 1
loadings[1,1] * ((evs[1]/var(climateSummerscaled$sum_meanT))^1/2)
loadings[2,1] * ((evs[1]/var(climateSummerscaled$sum_precip))^1/2)
loadings[3,1] * ((evs[1]/var(climateSummerscaled$sum_moisturedeficit))^1/2)
loadings[4,1] * ((evs[1]/var(climateSummerscaled$sum_PET))^1/2)
loadings[5,1] * ((evs[1]/var(climateSummerscaled$sum_GDD))^1/2)
loadings[6,1] * ((evs[1]/var(climateSummerscaled$fivedayrunning5C))^1/2)
loadings[7,1] * ((evs[1]/var(climateSummerscaled$fivedayrunning12C))^1/2)
loadings[8,1] * ((evs[1]/var(climateSummerscaled$GSLthreedayneg3C))^1/2)
loadings[9,1] * ((evs[1]/var(climateSummerscaled$iceoff_GL4))^1/2)
# Axis 2
loadings[1,2] * ((evs[2]/var(climateSummerscaled$sum_meanT))^1/2)
loadings[2,2] * ((evs[2]/var(climateSummerscaled$sum_precip))^1/2)
loadings[3,2] * ((evs[2]/var(climateSummerscaled$sum_moisturedeficit))^1/2)
loadings[4,2] * ((evs[2]/var(climateSummerscaled$sum_PET))^1/2)
loadings[5,2] * ((evs[2]/var(climateSummerscaled$sum_GDD))^1/2)
loadings[6,2] * ((evs[2]/var(climateSummerscaled$fivedayrunning5C))^1/2)
loadings[7,2] * ((evs[2]/var(climateSummerscaled$fivedayrunning12C))^1/2)
loadings[8,2] * ((evs[2]/var(climateSummerscaled$GSLthreedayneg3C))^1/2)
loadings[9,2] * ((evs[2]/var(climateSummerscaled$iceoff_GL4))^1/2)

# Subset the data into 10 year windows
w1 <- climateSummer[1:10,]
w2 <- climateSummer[2:11,]
w3 <- climateSummer[3:12,]
w4 <- climateSummer[4:13,]
w5 <- climateSummer[5:14,]
w6 <- climateSummer[6:15,]
w7 <- climateSummer[7:16,]
w8 <- climateSummer[8:17,]
w9 <- climateSummer[9:18,]
w10 <- climateSummer[10:19,]
w11 <- climateSummer[11:20,]
w12 <- climateSummer[12:21,]
w13 <- climateSummer[13:22,]
w14 <- climateSummer[14:23,]
w15 <- climateSummer[15:24,]
w16 <- climateSummer[16:25,]
w17 <- climateSummer[17:26,]
w18 <- climateSummer[18:27,]
w19 <- climateSummer[19:28,]
w20 <- climateSummer[20:29,]
w21 <- climateSummer[21:30,]
w22 <- climateSummer[22:31,]
w23 <- climateSummer[23:32,]
w24 <- climateSummer[24:33,]
w25 <- climateSummer[25:34,]
w26 <- climateSummer[26:35,]
w27 <- climateSummer[27:36,]
w28 <- climateSummer[28:37,]

# Execute PCA on each window and store results for later extraction
sumallPCA1 <- rda(na.exclude(w1[,2:ncol(w1)]), scale=T)
s1 <- summary(sumallPCA1)
sumallPCA2 <- rda(na.exclude(w2[,2:ncol(w2)]), scale=T)
s2 <- summary(sumallPCA2)
sumallPCA3 <- rda(na.exclude(w3[,2:ncol(w3)]), scale=T)
s3 <- summary(sumallPCA3)
sumallPCA4 <- rda(na.exclude(w4[,2:ncol(w4)]), scale=T)
s4 <- summary(sumallPCA4)
sumallPCA5 <- rda(na.exclude(w5[,2:ncol(w5)]), scale=T)
s5 <- summary(sumallPCA5)
sumallPCA6 <- rda(na.exclude(w6[,2:ncol(w6)]), scale=T)
s6 <- summary(sumallPCA6)
sumallPCA7 <- rda(na.exclude(w7[,2:ncol(w7)]), scale=T)
s7 <- summary(sumallPCA7)
sumallPCA8 <- rda(na.exclude(w8[,2:ncol(w8)]), scale=T)
s8 <- summary(sumallPCA8)
sumallPCA9 <- rda(na.exclude(w9[,2:ncol(w9)]), scale=T)
s9 <- summary(sumallPCA9)
sumallPCA10 <- rda(na.exclude(w10[,2:ncol(w10)]), scale=T)
s10 <- summary(sumallPCA10)
sumallPCA11 <- rda(na.exclude(w11[,2:ncol(w11)]), scale=T)
s11 <- summary(sumallPCA11)
sumallPCA12 <- rda(na.exclude(w12[,2:ncol(w12)]), scale=T)
s12 <- summary(sumallPCA12)
sumallPCA13 <- rda(na.exclude(w13[,2:ncol(w13)]), scale=T)
s13 <- summary(sumallPCA13)
sumallPCA14 <- rda(na.exclude(w14[,2:ncol(w14)]), scale=T)
s14 <- summary(sumallPCA14)
sumallPCA15 <- rda(na.exclude(w15[,2:ncol(w15)]), scale=T)
s15 <- summary(sumallPCA15)
sumallPCA16 <- rda(na.exclude(w16[,2:ncol(w16)]), scale=T)
s16 <- summary(sumallPCA16)
sumallPCA17 <- rda(na.exclude(w17[,2:ncol(w17)]), scale=T)
s17 <- summary(sumallPCA17)
sumallPCA18 <- rda(na.exclude(w18[,2:ncol(w18)]), scale=T)
s18 <- summary(sumallPCA18)
sumallPCA19 <- rda(na.exclude(w19[,2:ncol(w19)]), scale=T)
s19 <- summary(sumallPCA19)
sumallPCA20 <- rda(na.exclude(w20[,2:ncol(w20)]), scale=T)
s20 <- summary(sumallPCA20)
sumallPCA21 <- rda(na.exclude(w21[,2:ncol(w21)]), scale=T)
s21 <- summary(sumallPCA21)
sumallPCA22 <- rda(na.exclude(w22[,2:ncol(w22)]), scale=T)
s22 <- summary(sumallPCA22)
sumallPCA23 <- rda(na.exclude(w23[,2:ncol(w23)]), scale=T)
s23 <- summary(sumallPCA23)
sumallPCA24 <- rda(na.exclude(w24[,2:ncol(w24)]), scale=T)
s24 <- summary(sumallPCA24)
sumallPCA25 <- rda(na.exclude(w25[,2:ncol(w25)]), scale=T)
s25 <- summary(sumallPCA25)
sumallPCA26 <- rda(na.exclude(w26[,2:ncol(w26)]), scale=T)
s26 <- summary(sumallPCA26)
sumallPCA27 <- rda(na.exclude(w27[,2:ncol(w27)]), scale=T)
s27 <- summary(sumallPCA27)
sumallPCA28 <- rda(na.exclude(w28[,2:ncol(w28)]), scale=T)
s28 <- summary(sumallPCA28)



######################### Variable Score Analysis ##############################
# Make dataframes to store results of variable axes scores over the windows
# Axis 1
spec.scores.1 <- matrix(data = NA, nrow = 28, ncol = 10)
colnames(spec.scores.1) <- colnames(climateSummer)
spec.scores.1 <- as.data.frame(spec.scores.1)
colnames(spec.scores.1)[colnames(spec.scores.1)=="year"] <- "window"
spec.scores.1$window <- seq(1,28,1)

spec.scores.2 <- matrix(data = NA, nrow = 28, ncol = 10)
colnames(spec.scores.2) <- colnames(climateSummer)
spec.scores.2 <- as.data.frame(spec.scores.2)
colnames(spec.scores.2)[colnames(spec.scores.2)=="year"] <- "window"
spec.scores.2$window <- seq(1,28,1)

# Add the axis score to the dataframe
# Axis 1
for (i in 1:28) {
  spec.scores.1[i,2:10] <- get(paste("s",i,sep=""))[["species"]][,1]
}
# Axis 2
for (i in 1:28) {
  spec.scores.2[i,2:10] <- get(paste("s",i,sep=""))[["species"]][,2]
}

# Graph all of the ordinations into huge multipanel
par("mar")
par(mfrow = c(5,6),
    mar = c(1,1,1,1))
for (i in 1:28) {
  biplot(get(paste0("sumallPCA",i,sep="")))
}

# Now explore how the variable loadings change
# Make new long dataframe with both axes
a1 <- spec.scores.1[,2:10]
long1 <- melt(a1)
a2 <- spec.scores.2[,2:10]
long2 <- melt(a2)
long <- cbind(long1,long2)
long <- long[,c(1,2,4)]
colnames(long) <- c("Variable","PC1","PC2")

ggplot(long,aes(x = PC1, y = PC2)) + 
  geom_point(size = 3) +
  labs(y = "PC2",
       x = "PC1") +
  facet_wrap(~ Variable) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

# After looking at this, we can see that depending on the window of time, sometimes variables will still load strongly on an axis but they switch from negative or positive. Thus, let's take the absolute value, to focus on how the magnitude changes over time.
long$PC1 <- abs(long$PC1)
long$PC2 <- abs(long$PC2)

# Check the Axis magnitudes averaged across the windows and make table
mag1 <- as.data.frame(aggregate(long$PC1, by = list(long$Variable), FUN = mean))
mag1 <- mag1[order(mag1$x, decreasing = TRUE),]
colnames(mag1) <- c("Variable","PC1")
mag2 <- as.data.frame(aggregate(long$PC2, by = list(long$Variable), FUN = mean))
mag2 <- mag2[order(mag2$x, decreasing = TRUE),]
colnames(mag2) <- c("Variable","PC2")
#pdf(file = "Axis1WindowsMean.pdf", width = 3, height = 3)
#grid.table(mag1, rows = NULL)
#dev.off()
#pdf(file = "Axis2WindowsMean.pdf", width = 3, height = 3)
#grid.table(mag2, rows = NULL)
#dev.off()

# Compare to the entire dataset, make table
s.sp <- as.data.frame(s$species)
s.sp$PC1 <- abs(s.sp$PC1)
s.sp$PC2 <- abs(s.sp$PC2)
s.sp$Variable <- rownames(s.sp)
s.sp.1 <- as.data.frame(cbind(s.sp[,7],s.sp[,1]))
s.sp.2 <- as.data.frame(cbind(s.sp[,7],s.sp[,2]))
colnames(s.sp.1) <- c("Variable","PC1")
colnames(s.sp.2) <- c("Variable","PC2")
s.sp.1 <- s.sp.1[order(s.sp.1$PC1, decreasing = TRUE),]
s.sp.2 <- s.sp.2[order(s.sp.2$PC2, decreasing = TRUE),]
#pdf(file = "extended_summer/sensitivity_analyses/figs/Axis1AllData.pdf", width = 3, height = 3)
#grid.table(s.sp.1, rows = NULL)
#dev.off()
#pdf(file = "extended_summer/sensitivity_analyses/figs/Axis2AllData.pdf", width = 3, height = 3)
#grid.table(s.sp.2, rows = NULL)
#dev.off()

# Reorder by Axis 1 magnitude and rename variables
long$Variable = factor(long$Variable, levels=c("sum_meanT","sum_GDD","sum_PET","iceoff_GL4","fivedayrunning12C","fivedayrunning5C","sum_precip","sum_moisturedeficit","GSLthreedayneg3C"))
facet_names <-  c(`sum_meanT` = "Temp",`sum_GDD` = "GDD",`sum_PET` = "PET",`iceoff_GL4` = "Iceoff",`fivedayrunning12C` = "Days 12C",`fivedayrunning5C` = "Days 5C",`sum_precip` = "Precip",`sum_moisturedeficit` = "Moisture Deficit",`GSLthreedayneg3C` = "GSL")
ggplot(long,aes(x = PC1, y = PC2)) + 
  geom_point(size = 3, alpha = 0.25) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.75) +
  labs(y = "PC2",
       x = "PC1") +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

# Add window numbers to color for Caitlin
long$Window <- rep(1:28, times = 9)

# Graph with windows colored
ggplot(long,aes(x = PC1, y = PC2)) + 
  geom_point(pch = 16, size = 3, aes(colour = Window), alpha = 0.5) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.5, size = 0.25) +
  labs(y = "PC2",
       x = "PC1",
       title = "10-year windows (n = 28)") +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

min(long$PC1)
max(long$PC1)
min(long$PC2)
max(long$PC2)

min(longg$PC1)
max(longg$PC1)
min(longg$PC2)
max(longg$PC2)

# Figure 2 panel 1 (update with 2-color scale)
pdf("figs/Figure2a.pdf", width = 6, height = 4.5)
ggplot(long,aes(x = PC1, y = PC2)) + 
  geom_point(pch = 16, size = 2, aes(colour = Window)) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.5, size = 0.25) +
  labs(x = "PC1 Score of Variable",
       y = "PC2 Score of Variable",
       title = "a) 10-yr windows") +
  scale_x_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_gradient(low = "blue", high = "red") +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0, vjust = -0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()

# New Supplementary Figure - PC1 scores with Window on x-axis
pdf("figs/SuppFigure1a.pdf", width = 6, height = 4.5)
ggplot(long,aes(x = Window, y = PC1)) + 
  geom_point(pch = 16, size = 2, alpha = 0.75) +
  geom_line(alpha = 0.5, size = 0.25) +
  labs(x = "Window",
       y = "PC1 Score of Variable",
       title = "a) 10-yr windows") +
  scale_y_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0, vjust = -0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()



##################### Proportion Variation Explained Analysis ##################
# Make dataframes to store results
prop.expl <- matrix(data = NA, nrow = 28, ncol = 4)
colnames(prop.expl) <- c("Window","StartYear","PC1","PC2")
prop.expl <- as.data.frame(prop.expl)
prop.expl$Window <- seq(1,28,1)
prop.expl$StartYear <- seq(1982,2009,1)

# Add proportion explained data
for (i in 1:28) {
  prop.expl[i,3] <- get(paste("s",i,sep=""))[["cont"]][["importance"]][2,1]
  prop.expl[i,4] <- get(paste("s",i,sep=""))[["cont"]][["importance"]][2,2]
}

prop.expl$Comb <- prop.expl$PC1 + prop.expl$PC2

# Melt for ggplot
prop.expl.long <- melt(prop.expl,id.vars = c("Window","StartYear"),measure.vars = c("PC1","PC2"))
colnames(prop.expl.long) <- c("Window","StartYear","Axis","Prop")

# Plot
ggplot(prop.expl.long, aes(x = StartYear, y = Prop, shape = Axis)) + 
  geom_point(size = 3) +
  geom_line() +
  labs(y = "Proportion Variation Explained",
       x = "10-yr Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  theme_bw() +
  theme(legend.position = c(0.9,0.85),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1, units = "cm"),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))



############################ Year Score Analysis ###############################
# Make year score data frames for each window
for (i in 1:28) {
  assign(paste("yr",i,sep=""),as.data.frame(get(paste("s",i,sep=""))[["sites"]][,1:2]))
}

# Combine into one df and make year a column
year.scores <- rbind(yr1,yr2,yr3,yr4,yr5,yr6,yr7,yr8,yr9,yr10,yr11,yr12,yr13,yr14,yr15,yr16,yr17,yr18,yr19,yr20,yr21,yr22,yr23,yr24,yr25,yr26,yr27,yr28)
year.scores$Year <- rownames(year.scores)

# R altered the names of repeating years, so just extract the first 4 digits
year.scores$Year <- substr(year.scores$Year, 1, 4)
table(year.scores$Year)

# Subset to years with 10 windows of data
year.scores.10 <- subset(year.scores, Year == "1991"|Year == "1992"|Year == "1993"|Year == "1994"|Year == "1995"|Year == "1996"|Year == "1997"|Year == "1998"|Year == "1999"|Year == "2000"|Year == "2001"|Year == "2002"|Year == "2003"|Year == "2004"|Year == "2005"|Year == "2006"|Year == "2007"|Year == "2008"|Year == "2009")

# Take absolute values
year.scores.10$PC1 <- abs(year.scores.10$PC1)
year.scores.10$PC2 <- abs(year.scores.10$PC2)

# Look at which years are consistently "extended summers"
yraxis1 <- as.data.frame(aggregate(year.scores.10$PC1, by = list(year.scores.10$Year), FUN = mean))
yraxis1 <- yraxis1[order(yraxis1$x, decreasing = TRUE),]
yraxis1

# Graph
ggplot(year.scores.10,aes(x = PC1, y = PC2)) + 
  geom_point(size = 3, alpha = 0.25) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.75) +
  labs(y = "PC2",
       x = "PC1") +
  scale_x_continuous(labels = c(0,0.5,1,1.5,2)) +
  scale_y_continuous(labels = c(0,0.5,1,1.5,2)) +
  facet_wrap(~ Year) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black"),
        legend.key.size = unit(1.25, units = "cm"),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

# Update, for each year, get a moving window score
# e.g. for 1982 it will just be 1982, only 1 window. 1983 will be average of 2 windows. etc.
mean_mw_scores <- ddply(year.scores, "Year", summarise,
                        PC1 = mean(PC1),
                        PC2 = mean(PC2))
# How does it compare to the original scores
ggplot(data = NULL, aes(NWTclimate$sumallPC1, mean_mw_scores$PC1)) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = lm) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Original PC1",
       y = "Mean Window PC1")

ggplot(data = NULL, aes(NWTclimate$sumallPC2, mean_mw_scores$PC2)) +
  coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = lm) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Original PC2",
       y = "Mean Window PC2")

# PC1 general trend is the same. PC2 is very different! But not going to use PC2 anyway. The potential reason is that sometimes things driving PC2 like moisture may actually be relagated to PC3 in some windows! PC1 is always related to Temp and GDD.
# Save as file to test against saddle biomass
# write.csv(mean_mw_scores, "mean_window_scores.csv")


################################ Variable Removal ##############################
# Have the full PCA handy
# Full PCA
sumallPCA <- rda(na.exclude(climateSummer[,2:ncol(climateSummer)]), scale=T)
s <- summary(sumallPCA) 
s #PC1 (extended summer) explains 52.79%, PC2 (moisture gradient) explains 18.53%
par(mfrow = c(1,1),
    mar = c(4,4,4,4))
biplot(sumallPCA, main = "1982-ongoing")

# Remove Ice Off
climateSummer_noice <- climateSummer[,1:9]
sumallPCA_noice <- rda(na.exclude(climateSummer_noice[,2:ncol(climateSummer_noice)]), scale=T)
s_noice <- summary(sumallPCA_noice) 
s_noice # PC1 (extended summer) explains 53.81%, PC2 (moisture gradient) explains 20.58%
par(mfrow = c(1,1))
biplot(sumallPCA_noice, main = "No Ice Off")
par(mfrow = c(1,2),
    mar = c(3,2,3,1))
biplot(sumallPCA, main = "All Variables", ylim = c(-1.5,1.5), xlim = c(-2,2))
biplot(sumallPCA_noice, main = "No Ice Off", ylim = c(-1.5,1.5), xlim = c(-2,2))
# Basically nothing changes, except the axes actually explain a bit more variation

# Now try removing either the 12C or 5C variable (do both and see which is better)
# Remove 12C
climateSummer_noiceno12 <- climateSummer_noice[,-8]
sumallPCA_noiceno12 <- rda(na.exclude(climateSummer_noiceno12[,2:ncol(climateSummer_noiceno12)]), scale=T)
s_noiceno12 <- summary(sumallPCA_noiceno12) 
s_noiceno12 # PC1 (extended summer) explains 52.55%, PC2 (moisture gradient) explains 23.39%
par(mfrow = c(1,1))
biplot(sumallPCA_noiceno12, main = "No Ice Off or 12C")
par(mfrow = c(1,2),
    mar = c(3,2,3,1))
biplot(sumallPCA_noice, main = "No Ice Off", ylim = c(-1.5,1.5), xlim = c(-2,2))
biplot(sumallPCA_noiceno12, main = "No Ice Off or 12C", ylim = c(-1.5,1.5), xlim = c(-2,2))

# Remove 5C
climateSummer_noiceno5 <- climateSummer_noice[,-7]
sumallPCA_noiceno5 <- rda(na.exclude(climateSummer_noiceno5[,2:ncol(climateSummer_noiceno5)]), scale=T)
s_noiceno5 <- summary(sumallPCA_noiceno5) 
s_noiceno5 # PC1 (extended summer) explains 57.78%, PC2 (moisture gradient) explains 22.8%
par(mfrow = c(1,1))
biplot(sumallPCA_noiceno5, main = "No Ice Off or 5C")
par(mfrow = c(1,2),
    mar = c(3,2,3,1))
biplot(sumallPCA_noice, main = "No Ice Off", ylim = c(-1.5,1.5), xlim = c(-2,2))
biplot(sumallPCA_noiceno5, main = "No Ice Off or 5C", ylim = c(-1.5,1.5), xlim = c(-2,2))

# I'm not sure if we really need to remove variables. Some people have said the ratio of observations to variables should be 3:1 and we start with 4:1. The only reason to remove a variable is if others won't have it, such as Ice Off.

# Besides visual inspection, we could also try variable selection (though again this isn't really necessary because we only have 9 variables, some people may have hundreds or thousands)



############################# 15 year windows ##################################
# Repeat above analysis with 15 yr windows instead of 10 year windows
# Subset the data into 15 year windows
x1 <- climateSummer[1:15,]
x2 <- climateSummer[2:16,]
x3 <- climateSummer[3:17,]
x4 <- climateSummer[4:18,]
x5 <- climateSummer[5:19,]
x6 <- climateSummer[6:20,]
x7 <- climateSummer[7:21,]
x8 <- climateSummer[8:22,]
x9 <- climateSummer[9:23,]
x10 <- climateSummer[10:24,]
x11 <- climateSummer[11:25,]
x12 <- climateSummer[12:26,]
x13 <- climateSummer[13:27,]
x14 <- climateSummer[14:28,]
x15 <- climateSummer[15:29,]
x16 <- climateSummer[16:30,]
x17 <- climateSummer[17:31,]
x18 <- climateSummer[18:32,]
x19 <- climateSummer[19:33,]
x20 <- climateSummer[20:34,]
x21 <- climateSummer[21:35,]
x22 <- climateSummer[22:36,]
x23 <- climateSummer[23:37,]

# Execute PCA on each window and store results for later extraction
sumall.1 <- rda(na.exclude(x1[,2:ncol(x1)]), scale=T)
z1 <- summary(sumall.1)
sumall.2 <- rda(na.exclude(x2[,2:ncol(x2)]), scale=T)
z2 <- summary(sumall.2)
sumall.3 <- rda(na.exclude(x3[,2:ncol(x3)]), scale=T)
z3 <- summary(sumall.3)
sumall.4 <- rda(na.exclude(x4[,2:ncol(x4)]), scale=T)
z4 <- summary(sumall.4)
sumall.5 <- rda(na.exclude(x5[,2:ncol(x5)]), scale=T)
z5 <- summary(sumall.5)
sumall.6 <- rda(na.exclude(x6[,2:ncol(x6)]), scale=T)
z6 <- summary(sumall.6)
sumall.7 <- rda(na.exclude(x7[,2:ncol(x7)]), scale=T)
z7 <- summary(sumall.7)
sumall.8 <- rda(na.exclude(x8[,2:ncol(x8)]), scale=T)
z8 <- summary(sumall.8)
sumall.9 <- rda(na.exclude(x9[,2:ncol(x9)]), scale=T)
z9 <- summary(sumall.9)
sumall.10 <- rda(na.exclude(x10[,2:ncol(x10)]), scale=T)
z10 <- summary(sumall.10)
sumall.11 <- rda(na.exclude(x11[,2:ncol(x11)]), scale=T)
z11 <- summary(sumall.11)
sumall.12 <- rda(na.exclude(x12[,2:ncol(x12)]), scale=T)
z12 <- summary(sumall.12)
sumall.13 <- rda(na.exclude(x13[,2:ncol(x13)]), scale=T)
z13 <- summary(sumall.13)
sumall.14 <- rda(na.exclude(x14[,2:ncol(x14)]), scale=T)
z14 <- summary(sumall.14)
sumall.15 <- rda(na.exclude(x15[,2:ncol(x15)]), scale=T)
z15 <- summary(sumall.15)
sumall.16 <- rda(na.exclude(x16[,2:ncol(x16)]), scale=T)
z16 <- summary(sumall.16)
sumall.17 <- rda(na.exclude(x17[,2:ncol(x17)]), scale=T)
z17 <- summary(sumall.17)
sumall.18 <- rda(na.exclude(x18[,2:ncol(x18)]), scale=T)
z18 <- summary(sumall.18)
sumall.19 <- rda(na.exclude(x19[,2:ncol(x19)]), scale=T)
z19 <- summary(sumall.19)
sumall.20 <- rda(na.exclude(x20[,2:ncol(x20)]), scale=T)
z20 <- summary(sumall.20)
sumall.21 <- rda(na.exclude(x21[,2:ncol(x21)]), scale=T)
z21 <- summary(sumall.21)
sumall.22 <- rda(na.exclude(x22[,2:ncol(x22)]), scale=T)
z22 <- summary(sumall.22)
sumall.23 <- rda(na.exclude(x23[,2:ncol(x23)]), scale=T)
z23 <- summary(sumall.23)

#### Variable Score Analysis
# Make dataframes to store results of variable axes scores over the windows
# Axis 1
spec.scores.3 <- matrix(data = NA, nrow = 23, ncol = 10)
colnames(spec.scores.3) <- colnames(climateSummer)
spec.scores.3 <- as.data.frame(spec.scores.3)
colnames(spec.scores.3)[colnames(spec.scores.3)=="year"] <- "window"
spec.scores.3$window <- seq(1,23,1)

spec.scores.4 <- matrix(data = NA, nrow = 23, ncol = 10)
colnames(spec.scores.4) <- colnames(climateSummer)
spec.scores.4 <- as.data.frame(spec.scores.4)
colnames(spec.scores.4)[colnames(spec.scores.4)=="year"] <- "window"
spec.scores.4$window <- seq(1,23,1)

# Add the axis score to the dataframe
# Axis 1
for (i in 1:23) {
  spec.scores.3[i,2:10] <- get(paste("z",i,sep=""))[["species"]][,1]
}
# Axis 2
for (i in 1:23) {
  spec.scores.4[i,2:10] <- get(paste("z",i,sep=""))[["species"]][,2]
}

# Graph all of the ordinations into huge multipanel
par("mar")
par(mfrow = c(5,5),
    mar = c(1,1,1,1))
for (i in 1:23) {
  biplot(get(paste0("sumall.",i,sep="")))
}

# Now explore how the variable loadings change
# Make new long dataframe with both axes
ax1 <- spec.scores.3[,2:10]
longg1 <- melt(ax1)
ax2 <- spec.scores.4[,2:10]
longg2 <- melt(ax2)
longg <- cbind(longg1,longg2)
longg <- longg[,c(1,2,4)]
colnames(longg) <- c("Variable","PC1","PC2")

# Take the absolute value, to focus on how the magnitude changes over time.
longg$PC1 <- abs(longg$PC1)
longg$PC2 <- abs(longg$PC2)

# Check order of axis 1 variable strength
magg1 <- as.data.frame(aggregate(longg$PC1, by = list(longg$Variable), FUN = mean))
magg1 <- magg1[order(magg1$x, decreasing = TRUE),]
magg1 # Pretty similar to 10 year windows

# Keep same panel order to compare graphs
longg$Variable = factor(longg$Variable, levels=c("sum_meanT","sum_GDD","sum_PET","iceoff_GL4","fivedayrunning12C","fivedayrunning5C","sum_precip","sum_moisturedeficit","GSLthreedayneg3C"))

facet_names <-  c(`sum_meanT` = "Temp",`sum_GDD` = "GDD",`sum_PET` = "PET",`iceoff_GL4` = "Iceoff",`fivedayrunning12C` = "Days 12C",`fivedayrunning5C` = "Days 5C",`sum_precip` = "Precip",`sum_moisturedeficit` = "Moisture Deficit",`GSLthreedayneg3C` = "GSL")

# Add window numbers
longg$Window <- rep(1:23, times = 9)

# Plot
ggplot(longg, aes(x = PC1, y = PC2)) + 
  geom_point(pch = 16, size = 3, aes(colour = Window), alpha = 0.5) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.5, size = 0.25) +
  labs(y = "PC2",
       x = "PC1",
       title = "15-year windows (n = 23)") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

# Variation Explained
prop.expl <- matrix(data = NA, nrow = 23, ncol = 4)
colnames(prop.expl) <- c("Window","StartYear","PC1","PC2")
prop.expl <- as.data.frame(prop.expl)
prop.expl$Window <- seq(1,23,1)
prop.expl$StartYear <- seq(1982,2004,1)

# Add proportion explained data
for (i in 1:23) {
  prop.expl[i,3] <- get(paste("z",i,sep=""))[["cont"]][["importance"]][2,1]
  prop.expl[i,4] <- get(paste("z",i,sep=""))[["cont"]][["importance"]][2,2]
}

prop.expl$Comb <- prop.expl$PC1 + prop.expl$PC2

# Figure 2 panel 1 (update with 2-color scale)
pdf("figs/Figure2b.pdf", width = 6, height = 4.5)
ggplot(longg, aes(x = PC1, y = PC2)) + 
  geom_point(pch = 16, size = 2, aes(colour = Window)) +
  geom_path(arrow = arrow(angle = 20, length = unit(0.2,"cm")), alpha = 0.5, size = 0.25) +
  labs(x = "PC1 Score of Variable",
       y = "PC2 Score of Variable",
       title = "b) 15-yr windows") +
  scale_x_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_gradient(low = "blue", high = "red") +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0, vjust = -0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()

# New Supplementary Figure - PC1 scores with Window on x-axis
pdf("figs/SuppFigure1b.pdf", width = 6, height = 4.5)
ggplot(longg, aes(x = Window, y = PC1)) + 
  geom_point(pch = 16, size = 2, alpha = 0.75) +
  geom_line(alpha = 0.5, size = 0.25) +
  labs(x = "Window",
       y = "PC1 Score of Variable",
       title = "b) 15-yr windows") +
  scale_y_continuous(limits = c(0,1.1),
                     breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     labels = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ Variable, labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0, vjust = -0.5),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()



########################### Check D1 Correlations ##############################
# Above analyses show GSL doesn't always load strongly on PC1. We want to see if it is related to winter precipitation. So, look at correlations with D1 precip data
d1 <- read.csv("../../climate_d1_c1/output_data/d1_dailyppt_infilled_1952-2018.csv")
# Now need to summarize for each year, sum of precip. So year 1982 will be Nov 1981 to May 1982.
d1 <- subset(d1, year > 1980)
climateSummer$d1wp <- NA
for (i in 1982:2018) {
  sub1 <- subset(d1, year == i - 1 & month > 10)
  sum1 <- sum(sub1$precipitation.mm)
  sub2 <- subset(d1, year == i & month < 6)
  sum2 <- sum(sub2$precipitation.mm)
  climateSummer$d1wp[i-1981] <- sum(sum1, sum2)
}

# Okay, now we have D1 winter precip as a column with the other climate summer variables
# Make the 28 10-year windows
w1 <- climateSummer[1:10,]
w2 <- climateSummer[2:11,]
w3 <- climateSummer[3:12,]
w4 <- climateSummer[4:13,]
w5 <- climateSummer[5:14,]
w6 <- climateSummer[6:15,]
w7 <- climateSummer[7:16,]
w8 <- climateSummer[8:17,]
w9 <- climateSummer[9:18,]
w10 <- climateSummer[10:19,]
w11 <- climateSummer[11:20,]
w12 <- climateSummer[12:21,]
w13 <- climateSummer[13:22,]
w14 <- climateSummer[14:23,]
w15 <- climateSummer[15:24,]
w16 <- climateSummer[16:25,]
w17 <- climateSummer[17:26,]
w18 <- climateSummer[18:27,]
w19 <- climateSummer[19:28,]
w20 <- climateSummer[20:29,]
w21 <- climateSummer[21:30,]
w22 <- climateSummer[22:31,]
w23 <- climateSummer[23:32,]
w24 <- climateSummer[24:33,]
w25 <- climateSummer[25:34,]
w26 <- climateSummer[26:35,]
w27 <- climateSummer[27:36,]
w28 <- climateSummer[28:37,]

# Correlations
cor.df <- as.data.frame(matrix(NA, 28, 6))
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["d1wp"]],get(paste("w",i,sep=""))[["GSLthreedayneg3C"]])
  cor.df[i,1] <- c$estimate
  cor.df[i,2] <- c$p.value
}
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["d1wp"]],get(paste("w",i,sep=""))[["sum_meanT"]])
  cor.df[i,3] <- c$estimate
  cor.df[i,4] <- c$p.value
}
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["GSLthreedayneg3C"]],get(paste("w",i,sep=""))[["sum_meanT"]])
  cor.df[i,5] <- c$estimate
  cor.df[i,6] <- c$p.value
}
colnames(cor.df) <- c("P.GSL.R", "P.GSL.P","P.T.R","P.T.P","GSL.T.R","GSL.T.P")

# Make significance dataframe
sig.df <- as.data.frame(matrix(NA, 28, 3))
colnames(sig.df) <- c("P.GSL.Sig","P.T.Sig","GSL.T.Sig")
for (i in 1:nrow(cor.df)) {
  if (cor.df$P.GSL.P[i] < 0.05) {
    sig.df$P.GSL.Sig[i] <- "Significant"
  } else {
    sig.df$P.GSL.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$P.T.P[i] < 0.05) {
    sig.df$P.T.Sig[i] <- "Significant"
  } else {
    sig.df$P.T.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$GSL.T.P[i] < 0.05) {
    sig.df$GSL.T.Sig[i] <- "Significant"
  } else {
    sig.df$GSL.T.Sig[i] <- "Not Significant"
  }
}



# Melt for plotting
cor.df.long <- melt(cor.df)
sig.df.long <- melt(sig.df, id.vars = NULL)

# Define r or p value
cor.df.long$stat[c(29:56,85:112,141:168)] <- "p"
cor.df.long$stat[c(1:28,57:84,113:140)] <- "r"

# Get just the r's but can color by significance
cor.df.long.r <- subset(cor.df.long, stat == "r")

# Add significant column
cor.df.long.r$sig <- sig.df.long$value

# Add the windows for the x axis
cor.df.long.r$window <- rep(1:28,time = 3)
cor.df.long.r$windowstart <- rep(1982:2009, time = 3)

# Plot
cor.df.long.r$variable <- factor(cor.df.long.r$variable, levels = c("GSL.T.R","P.GSL.R","P.T.R"))
ggplot(cor.df.long.r,aes(x=windowstart,y=value,colour=variable,shape=sig,group=variable)) + 
  geom_hline(yintercept = 0) +
  geom_point(size = 3) +
  geom_line() +
  labs(y = "Correlation Coefficient",
       x = "10-yr Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  scale_colour_discrete(labels = c("GSL, Temp","Snow, GSL", "Snow, Temp")) +
  scale_shape_discrete(labels = c("p > 0.05","p < 0.05")) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(1, units = "cm"),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))

# Caitlin says that data artifacts may be driving trends seen in the above graph. To check, let's do the same correlations as above except with D1 temp instead of saddle temp.
d1t <- read.csv("../../climate_d1_c1/output_data/d1_dailytemp_infilled_1952-2018.csv")
# Now need to summarize for each year, mean of daily mean temp for Jun - Aug
d1t <- subset(d1t, year > 1981)
climateSummer$d1temp <- NA
for (i in 1982:2018) {
  sub1 <- subset(d1t, year == i & month==6)
  sub2 <- subset(d1t, year == i & month==7)
  sub3 <- subset(d1t, year == i & month==8)
  summer <- rbind(sub1,sub2,sub3)
  climateSummer$d1temp[i-1981] <- mean(summer$Tmean.C)
}

# gsl is number of days bounded by 3 consec days tmin < -3C 
# emily split search windows on july 15 (jday = 196): latest date in spring (jday < 196) that 3day consec tmin -3 & first date in fall (jday > 196) that 3day consec tmin -3
## expect for d1 will be shorter than sdl temp, bc higher elevation (cooler), but also sdl is electronic logger, which typically records warmer temps than chart recorders (but also colder tmins than chart)
d1_gsl <- d1t %>%
  # logical: is tmin < -3?
  mutate(neg3 = Tmin.C < -3,
         # tally number of days in 3-day window that tmin < -3
         consecneg3 = (lag(neg3, 1) + lag(neg3, 2) + neg3)) %>%
  select(year, doy, Tmin.C, neg3, consecneg3) %>%
  group_by(year) %>%
  # ID latest spring date that had consec 3-day tmin < -3
  mutate(sprdoy = max(doy[doy < 196 & consecneg3 == 3], na.rm = T),
         # ID earliest fall date that had consec 3-day tmin < -3
         faldoy = min(doy[doy > 196 & consecneg3 == 3], na.rm = T),
         # diff fal from spring date to get GSL
         d1neg3GSL = faldoy-sprdoy) %>%
  ungroup() %>%
  # winnow df down to yearly GSL values
  select(year, sprdoy, faldoy, d1neg3GSL) %>%
  distinct()

# join d1 GSL to climateSummer
climateSummer <- merge(climateSummer, d1_gsl[c("year", "d1neg3GSL")]) 

w1 <- climateSummer[1:10,]
w2 <- climateSummer[2:11,]
w3 <- climateSummer[3:12,]
w4 <- climateSummer[4:13,]
w5 <- climateSummer[5:14,]
w6 <- climateSummer[6:15,]
w7 <- climateSummer[7:16,]
w8 <- climateSummer[8:17,]
w9 <- climateSummer[9:18,]
w10 <- climateSummer[10:19,]
w11 <- climateSummer[11:20,]
w12 <- climateSummer[12:21,]
w13 <- climateSummer[13:22,]
w14 <- climateSummer[14:23,]
w15 <- climateSummer[15:24,]
w16 <- climateSummer[16:25,]
w17 <- climateSummer[17:26,]
w18 <- climateSummer[18:27,]
w19 <- climateSummer[19:28,]
w20 <- climateSummer[20:29,]
w21 <- climateSummer[21:30,]
w22 <- climateSummer[22:31,]
w23 <- climateSummer[23:32,]
w24 <- climateSummer[24:33,]
w25 <- climateSummer[25:34,]
w26 <- climateSummer[26:35,]
w27 <- climateSummer[27:36,]
w28 <- climateSummer[28:37,]

# Correlations
cor.df <- as.data.frame(matrix(NA, 28, 6))
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["d1wp"]],get(paste("w",i,sep=""))[["d1neg3GSL"]])
  cor.df[i,1] <- c$estimate
  cor.df[i,2] <- c$p.value
}
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["d1wp"]],get(paste("w",i,sep=""))[["d1temp"]])
  cor.df[i,3] <- c$estimate
  cor.df[i,4] <- c$p.value
}
for (i in 1:28) {
  c <- cor.test(get(paste("w",i,sep=""))[["d1neg3GSL"]],get(paste("w",i,sep=""))[["d1temp"]])
  cor.df[i,5] <- c$estimate
  cor.df[i,6] <- c$p.value
}
colnames(cor.df) <- c("P.GSL.R", "P.GSL.P","P.T.R","P.T.P","GSL.T.R","GSL.T.P")

# Make significance dataframe
sig.df <- as.data.frame(matrix(NA, 28, 3))
colnames(sig.df) <- c("P.GSL.Sig","P.T.Sig","GSL.T.Sig")
for (i in 1:nrow(cor.df)) {
  if (cor.df$P.GSL.P[i] < 0.05) {
    sig.df$P.GSL.Sig[i] <- "Significant"
  } else {
    sig.df$P.GSL.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$P.T.P[i] < 0.05) {
    sig.df$P.T.Sig[i] <- "Significant"
  } else {
    sig.df$P.T.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$GSL.T.P[i] < 0.05) {
    sig.df$GSL.T.Sig[i] <- "Significant"
  } else {
    sig.df$GSL.T.Sig[i] <- "Not Significant"
  }
}



# Melt for plotting
cor.df.long <- melt(cor.df)
sig.df.long <- melt(sig.df, id.vars = NULL)

# Define r or p value
cor.df.long$stat[c(29:56,85:112,141:168)] <- "p"
cor.df.long$stat[c(1:28,57:84,113:140)] <- "r"

# Get just the r's but can color by significance
cor.df.long.r <- subset(cor.df.long, stat == "r")

# Add significant column
cor.df.long.r$sig <- sig.df.long$value

# Add the windows for the x axis
cor.df.long.r$window <- rep(1:28,time = 3)
cor.df.long.r$windowstart <- rep(1982:2009, time = 3)

# Plot
cor.df.long.r$variable <- factor(cor.df.long.r$variable, levels = c("GSL.T.R","P.GSL.R","P.T.R"))
ggplot(cor.df.long.r,aes(x=windowstart,y=value,colour=variable,shape=sig,group=variable)) + 
  geom_hline(yintercept = 0) +
  geom_point(size = 3) +
  geom_line() +
  labs(y = "Correlation Coefficient",
       x = "10-yr Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  scale_colour_discrete(labels = c("GSL, Temp","Snow, GSL", "Snow, Temp")) +
  scale_shape_discrete(labels = c("p > 0.05","p < 0.05")) +
  labs(title = "D1: Summer Tmean, Oct-May ppt, and GSL correlations",
       subtitle = "Last several GSL/Temp points potentially influenced by artificial drops in D1 chart T") +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(1, units = "cm"),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18),
        plot.subtitle = element_text(size = 10))

# Plot for Supp Fig 3
sf1 <- subset(cor.df.long.r, variable == "GSL.T.R")
pdf("figs/SuppFig3.pdf", width = 6.7, height = 4.55)
ggplot(data = sf1) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = value, shape = sig),colour = "blue",size = 4) +
  geom_line(aes(x = windowstart, y = value),colour = "blue",size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("p > 0.05","p < 0.05")) +
  scale_y_continuous(breaks = c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
                     limits = c(-0.75, 1)) +
  theme(legend.position = "none",
        axis.title.x = element_text(face="bold", size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(face="bold",size=18))
dev.off()

# Sort of similar to SDL trends.. red line hump less pronounced in early windows and spiked drop in latter windows (but could be from data/instrument QA issues .. there were 3 summers with artificial drops in temp)
# snow v temp and snow v GSL are more consistent
# ggsave("extended_summer/sensitivity_analyses/figs/CorrelationsD1temp.pdf", width = 6, height = 4, units = "in")
# ctw: difference in 80s windows btwn SDL and D1, in corr values and signif levels is problematic to me.. but see what Katie thinks
# may not need to be concerned about exact correlations, more than they flux through time depending on time window used
# also, ctw is now noticing most of the point in SDL and all in D1 are NOT signif, so I care less now..
# on the other hand.. it's a cor test on n=10 so don't expect small pvals given weather does bop interannually

# consider D1 GSL and Tmean vs. @ SDL. I think logger projection started at SDL 1999 and backwards (electronic starts in 2000). Here, D1 data are all from D1 chart, no logger data
select(climateSummer, year, sum_meanT, d1temp, GSLthreedayneg3C, d1neg3GSL) %>%
  gather(met, val, sum_meanT:ncol(.)) %>%
  mutate(site = ifelse(grepl("d1", met), "D1", "SDL"),
         met2 = ifelse(grepl("GSL", met), "GSL", "summer meanT")) %>%
  ggplot(aes(year, val, col = site)) +
  geom_point(alpha = 0.6) +
  geom_smooth(aes(fill = site)) +
  facet_wrap(~met2, scales = "free_y")



################### New Correlation Analyses for Paper ########################
# Temp and GDD. Temp and GSL. To exemplify 2 strong axis 1 correlations, and 1 strong 1 weak
# Just saddle summer data used in PCA
cor.df <- as.data.frame(matrix(NA, 28, 6))
colnames(cor.df) <- c("T.GDD.R","T.GSL.R","T.GDD.P","T.GSL.P","T.GDD.Sig","T.GSL.Sig")
for (i in 1:nrow(cor.df)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_GDD"]])
  cor.df[i,1] <- c$estimate
  cor.df[i,3] <- c$p.value
}
for (i in 1:nrow(cor.df)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["GSLthreedayneg3C"]])
  cor.df[i,2] <- c$estimate
  cor.df[i,4] <- c$p.value
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$T.GDD.P[i] < 0.05) {
    cor.df$T.GDD.Sig[i] <- "Significant"
  } else {
    cor.df$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$T.GSL.P[i] < 0.05) {
    cor.df$T.GSL.Sig[i] <- "Significant"
  } else {
    cor.df$T.GSL.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
cor.df$windowstart <- rep(1982:2009, time = 1)

# Plot
pdf("figs/Correlations10yr.pdf", width = 4, height = 3)
ggplot(data = cor.df) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 1990, linetype = "dotted") +
  geom_vline(xintercept = 1997, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GDD.R, shape = T.GDD.Sig),colour = "red",size = 2) +
  geom_line(aes(x = windowstart, y = T.GDD.R),colour = "red",size = 0.5) +
  geom_point(aes(x = windowstart, y = T.GSL.R, shape = T.GSL.Sig),colour = "blue",size = 2) +
  geom_line(aes(x = windowstart, y = T.GSL.R),colour = "blue",size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  scale_colour_manual(values = c("blue","red"),
                      labels = c("Temp, GDD","Temp, GSL")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("p > 0.05","p < 0.05")) +
  theme_bw() +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1),
                     limits = c(-0.25, 1)) +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()

# Breakpoint regression of window start year and correlation coefficient
bp.df.10 <- as.data.frame(cor.df$T.GSL.R)
CE.Normal.Mean(bp.df.10, penalty = "AIC", distyp = 1)
CE.Normal.Mean(bp.df.10, penalty = "BIC", distyp = 1)
CE.Normal.Mean(bp.df.10, penalty = "AIC", distyp = 2)
CE.Normal.Mean(bp.df.10, penalty = "BIC", distyp = 2)
# Says there's one break point at 13 (window start 1994)
CE.Normal.MeanVar(bp.df.10, penalty = "AIC", distyp = 1)
CE.Normal.MeanVar(bp.df.10, penalty = "BIC", distyp = 1)
CE.Normal.MeanVar(bp.df.10, penalty = "AIC", distyp = 2)
CE.Normal.MeanVar(bp.df.10, penalty = "BIC", distyp = 2)
# Says there are two break points at 13, 22 (window start 1994, 2003)

# Segmented
m <- lm(T.GSL.R ~ windowstart, data = cor.df)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1994)))
summary(my.seg1)
my.seg2 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994)))
summary(my.seg2)
my.seg3 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994, 1998)))
summary(my.seg3)
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 2 breakpoints - 1990, 1997


### 15 yr windows
cor.df.15 <- as.data.frame(matrix(NA, 23, 6))
colnames(cor.df.15) <- c("T.GDD.R","T.GSL.R","T.GDD.P","T.GSL.P","T.GDD.Sig","T.GSL.Sig")
for (i in 1:nrow(cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_GDD"]])
  cor.df.15[i,1] <- c$estimate
  cor.df.15[i,3] <- c$p.value
}
for (i in 1:nrow(cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["GSLthreedayneg3C"]])
  cor.df.15[i,2] <- c$estimate
  cor.df.15[i,4] <- c$p.value
}
for (i in 1:nrow(cor.df.15)) {
  if (cor.df.15$T.GDD.P[i] < 0.05) {
    cor.df.15$T.GDD.Sig[i] <- "Significant"
  } else {
    cor.df.15$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df.15)) {
  if (cor.df.15$T.GSL.P[i] < 0.05) {
    cor.df.15$T.GSL.Sig[i] <- "Significant"
  } else {
    cor.df.15$T.GSL.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
cor.df.15$windowstart <- rep(1982:2004, time = 1)

# Plot
pdf("figs/Correlations15yr.pdf", width = 4, height = 3)
ggplot(data = cor.df.15) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 1992, linetype = "dotted") +
  geom_vline(xintercept = 1994, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GDD.R, shape = T.GDD.Sig),colour = "red",size = 2) +
  geom_line(aes(x = windowstart, y = T.GDD.R),colour = "red",size = 0.5) +
  geom_point(aes(x = windowstart, y = T.GSL.R, shape = T.GSL.Sig),colour = "blue",size = 2) +
  geom_line(aes(x = windowstart, y = T.GSL.R),colour = "blue",size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1989, 1996, 2003)) +
  scale_colour_manual(values = c("blue","red"),
                      labels = c("Temp, GDD","Temp, GSL")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("p > 0.05","p < 0.05")) +
  theme_bw() +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1),
                     limits = c(-0.25, 1)) +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()


# Breakpoint regression of window start year and correlation coefficient
bp.df.15 <- as.data.frame(cor.df.15$T.GSL.R)
CE.Normal.Mean(bp.df.15, penalty = "AIC", distyp = 1)
CE.Normal.Mean(bp.df.15, penalty = "BIC", distyp = 1)
CE.Normal.Mean(bp.df.15, penalty = "AIC", distyp = 2)
CE.Normal.Mean(bp.df.15, penalty = "BIC", distyp = 2)
# Says there's one break point at 13 (window start 1994)
CE.Normal.MeanVar(bp.df.15, penalty = "AIC", distyp = 1)
CE.Normal.MeanVar(bp.df.15, penalty = "BIC", distyp = 1)
CE.Normal.MeanVar(bp.df.15, penalty = "AIC", distyp = 2)
CE.Normal.MeanVar(bp.df.15, penalty = "BIC", distyp = 2)
# Says there are two break points at 12, 18 (window start 1994, 2003)

# Segmented (I think this is best)
m <- lm(T.GSL.R ~ windowstart, data = cor.df)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1994)))
summary(my.seg1)
my.seg2 <- segmented(m,
                    seg.Z = ~ windowstart, 
                    psi = list(windowstart = c(1992, 1994)))
summary(my.seg2)
my.seg3 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994, 1998)))
summary(my.seg3)
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 2 breakpoints - 1992, 1994

# Try strucchange
library(strucchange)
m1 <- breakpoints(T.GSL.R ~ windowstart, data = cor.df, breaks = 1, h = 3)
summary(m1)
m2 <- breakpoints(T.GSL.R ~ windowstart, data = cor.df, breaks = 2, h = 3)
summary(m2)

# Use these data for little insets in Figure 6
ggplot(data = cor.df) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GDD.R), size = 0.5) +
  labs(x = "Time",
       y = "r") +
  theme_classic() +
  scale_y_continuous(breaks = c(-1,0,1),
                     limits = c(-1, 1)) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

ggplot(data = cor.df) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GSL.R), size = 0.5) +
  labs(x = "Time",
       y = "r") +
  theme_classic() +
  scale_y_continuous(breaks = c(-1,0,1),
                     limits = c(-1, 1)) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

# Make little raw scatter plots to include on Figure 4
pdf("figs/mini_TGDD.pdf", width = 1.5, height = 1.5)
ggplot(data = climateSummer, aes(sum_meanT, sum_GDD)) + 
  geom_point(size = 1) +
  labs(x = "Temp",
       y = "GDD") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12))
dev.off()

pdf("figs/mini_TGSL.pdf", width = 1.5, height = 1.5)
ggplot(data = climateSummer, aes(sum_meanT, GSLthreedayneg3C)) + 
  geom_point(size = 1) +
  labs(x = "Temp",
       y = "GSL") +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 12))
dev.off()

### Emily wanted to check if randomly selecting 10 years from the whole dataset could yield 0 correlation between GSL and Temp or if this was really a moving window thing.
yrTgsl <- select(climateSummer, year, sum_meanT, GSLthreedayneg3C)
names(yrTgsl) <- c("Year","Temp","GSL")
yrTgsl$Decade <- "four"
for (i in 1:nrow(yrTgsl)) {
  if(yrTgsl$Year[i] < 2010) {
    yrTgsl$Decade[i] <- "three"
  }
}
for (i in 1:nrow(yrTgsl)) {
  if(yrTgsl$Year[i] < 2000) {
    yrTgsl$Decade[i] <- "two"
  }
}
for (i in 1:nrow(yrTgsl)) {
  if(yrTgsl$Year[i] < 1990) {
    yrTgsl$Decade[i] <- "one"
  }
}
yrTgsl$Decade <- as.factor(yrTgsl$Decade)
yrTgsl <- yrTgsl %>% group_by(Decade)
results_df <- as.data.frame(matrix(NA, 500, 3))
colnames(results_df) <- c("Iteration","R","P")
set.seed(500)
results_df$Iteration <- seq(1:500)
for (i in 1:nrow(results_df)) {
  d <- sample_n(yrTgsl, 3)
  c <- cor.test(d$Temp, d$GSL)
  results_df$R[i] <- c$estimate
  results_df$P[i] <- c$p.value
}
min(results_df$R) # 0.4. 
# Thus, when using 12 random samples from across the whole dataset (3 in each decade), r never goes down to zero, as seen in some of the moving windows



#### _ All variables vs. Temp ####
all.cor.df.10 <- as.data.frame(matrix(NA, 28, 24))
colnames(all.cor.df.10) <- c("T.GDD.R","T.GSL.R","T.Precip.R","T.MD.R","T.PET.R","T.five.R","T.twelve.R","T.iceoff.R","T.GDD.P","T.GSL.P","T.Precip.P","T.MD.P","T.PET.P","T.five.P","T.twelve.P","T.iceoff.P","T.GDD.Sig","T.GSL.Sig","T.Precip.Sig","T.MD.Sig","T.PET.Sig","T.five.Sig","T.twelve.Sig","T.iceoff.Sig")
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_GDD"]])
  all.cor.df.10[i,1] <- c$estimate
  all.cor.df.10[i,9] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["GSLthreedayneg3C"]])
  all.cor.df.10[i,2] <- c$estimate
  all.cor.df.10[i,10] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_precip"]])
  all.cor.df.10[i,3] <- c$estimate
  all.cor.df.10[i,11] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_moisturedeficit"]])
  all.cor.df.10[i,4] <- c$estimate
  all.cor.df.10[i,12] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_PET"]])
  all.cor.df.10[i,5] <- c$estimate
  all.cor.df.10[i,13] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["fivedayrunning5C"]])
  all.cor.df.10[i,6] <- c$estimate
  all.cor.df.10[i,14] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["fivedayrunning12C"]])
  all.cor.df.10[i,7] <- c$estimate
  all.cor.df.10[i,15] <- c$p.value
}
for (i in 1:nrow(all.cor.df.10)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["iceoff_GL4"]])
  all.cor.df.10[i,8] <- c$estimate
  all.cor.df.10[i,16] <- c$p.value
}


for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.GDD.P[i] < 0.05) {
    all.cor.df.10$T.GDD.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.GSL.P[i] < 0.05) {
    all.cor.df.10$T.GSL.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.GSL.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.Precip.P[i] < 0.05) {
    all.cor.df.10$T.Precip.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.Precip.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.MD.P[i] < 0.05) {
    all.cor.df.10$T.MD.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.MD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.PET.P[i] < 0.05) {
    all.cor.df.10$T.PET.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.PET.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.five.P[i] < 0.05) {
    all.cor.df.10$T.five.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.five.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.twelve.P[i] < 0.05) {
    all.cor.df.10$T.twelve.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.twelve.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.10)) {
  if (all.cor.df.10$T.iceoff.P[i] < 0.05) {
    all.cor.df.10$T.iceoff.Sig[i] <- "Significant"
  } else {
    all.cor.df.10$T.iceoff.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
all.cor.df.10$windowstart <- rep(1982:2009, time = 1)

# Organize dataframe and melt
all.cor.df.10 <- select(all.cor.df.10, -c(T.GDD.R, T.GSL.R,
                                          T.GDD.P, T.GSL.P, T.Precip.P, T.MD.P, T.PET.P,
                                          T.five.P, T.twelve.P, T.iceoff.P,
                                          T.GDD.Sig, T.GSL.Sig))
all.cor.long.10 <- melt(all.cor.df.10,
                        measure.vars = c("T.Precip.R", "T.MD.R", "T.PET.R", "T.five.R",
                                         "T.twelve.R", "T.iceoff.R",
                                         "T.Precip.Sig", "T.MD.Sig", "T.PET.Sig", "T.five.Sig",
                                         "T.twelve.Sig", "T.iceoff.Sig"))
cols <- all.cor.long.10[169:nrow(all.cor.long.10),]
all.cor.long.10 <- all.cor.long.10[1:168,]
all.cor.long.10.2 <- cbind(all.cor.long.10, cols)
all.cor.long.10.2 <- all.cor.long.10.2[,-(c(4,5))]
names(all.cor.long.10.2)[4] <- "sig"
all.cor.long.10.2$windowlength <- "a) 10-yr windows"

# 15 yr windows
all.cor.df.15 <- as.data.frame(matrix(NA, 23, 24))
colnames(all.cor.df.15) <- c("T.GDD.R","T.GSL.R","T.Precip.R","T.MD.R","T.PET.R","T.five.R","T.twelve.R","T.iceoff.R","T.GDD.P","T.GSL.P","T.Precip.P","T.MD.P","T.PET.P","T.five.P","T.twelve.P","T.iceoff.P","T.GDD.Sig","T.GSL.Sig","T.Precip.Sig","T.MD.Sig","T.PET.Sig","T.five.Sig","T.twelve.Sig","T.iceoff.Sig")
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_GDD"]])
  all.cor.df.15[i,1] <- c$estimate
  all.cor.df.15[i,9] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["GSLthreedayneg3C"]])
  all.cor.df.15[i,2] <- c$estimate
  all.cor.df.15[i,10] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_precip"]])
  all.cor.df.15[i,3] <- c$estimate
  all.cor.df.15[i,11] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_moisturedeficit"]])
  all.cor.df.15[i,4] <- c$estimate
  all.cor.df.15[i,12] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_PET"]])
  all.cor.df.15[i,5] <- c$estimate
  all.cor.df.15[i,13] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["fivedayrunning5C"]])
  all.cor.df.15[i,6] <- c$estimate
  all.cor.df.15[i,14] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["fivedayrunning12C"]])
  all.cor.df.15[i,7] <- c$estimate
  all.cor.df.15[i,15] <- c$p.value
}
for (i in 1:nrow(all.cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["iceoff_GL4"]])
  all.cor.df.15[i,8] <- c$estimate
  all.cor.df.15[i,16] <- c$p.value
}


for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.GDD.P[i] < 0.05) {
    all.cor.df.15$T.GDD.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.GSL.P[i] < 0.05) {
    all.cor.df.15$T.GSL.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.GSL.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.Precip.P[i] < 0.05) {
    all.cor.df.15$T.Precip.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.Precip.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.MD.P[i] < 0.05) {
    all.cor.df.15$T.MD.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.MD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.PET.P[i] < 0.05) {
    all.cor.df.15$T.PET.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.PET.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.five.P[i] < 0.05) {
    all.cor.df.15$T.five.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.five.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.twelve.P[i] < 0.05) {
    all.cor.df.15$T.twelve.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.twelve.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(all.cor.df.15)) {
  if (all.cor.df.15$T.iceoff.P[i] < 0.05) {
    all.cor.df.15$T.iceoff.Sig[i] <- "Significant"
  } else {
    all.cor.df.15$T.iceoff.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
all.cor.df.15$windowstart <- rep(1982:2004, time = 1)

# Organize dataframe and melt
all.cor.df.15 <- select(all.cor.df.15, -c(T.GDD.R, T.GSL.R,
                                          T.GDD.P, T.GSL.P, T.Precip.P, T.MD.P, T.PET.P,
                                          T.five.P, T.twelve.P, T.iceoff.P,
                                          T.GDD.Sig, T.GSL.Sig))
all.cor.long.15 <- melt(all.cor.df.15,
                        measure.vars = c("T.Precip.R", "T.MD.R", "T.PET.R", "T.five.R",
                                         "T.twelve.R", "T.iceoff.R",
                                         "T.Precip.Sig", "T.MD.Sig", "T.PET.Sig", "T.five.Sig",
                                         "T.twelve.Sig", "T.iceoff.Sig"))
cols.15 <- all.cor.long.15[139:nrow(all.cor.long.15),]
all.cor.long.15 <- all.cor.long.15[1:138,]
all.cor.long.15.2 <- cbind(all.cor.long.15, cols.15)
all.cor.long.15.2 <- all.cor.long.15.2[,-(c(4,5))]
names(all.cor.long.15.2)[4] <- "sig"
all.cor.long.15.2$windowlength <- "b) 15-yr windows"

# Combine and plot
all.cor.multi <- rbind(all.cor.long.10.2, all.cor.long.15.2)
all.cor.multi$value <- as.numeric(all.cor.multi$value)
all.cor.multi <- droplevels(all.cor.multi)
levels(all.cor.multi$variable) <- c("Precip","MD","PET","Days5C","Days12C","Iceoff")

pdf("figs/SuppFigure2.pdf", width = 7, height = 3)
ggplot(data = all.cor.multi, aes(windowstart, value, colour = variable)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_point(size = 2, aes(shape = sig)) +
  geom_line(size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  facet_wrap(~ windowlength, ncol = 2) +
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.spacing.y = unit(-0.1, "cm"),
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()



#### _Check Spearman correlations ####
cor.df <- as.data.frame(matrix(NA, 28, 6))
colnames(cor.df) <- c("T.GDD.R","T.GSL.R","T.GDD.P","T.GSL.P","T.GDD.Sig","T.GSL.Sig")
for (i in 1:nrow(cor.df)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["sum_GDD"]], method = "spearman")
  cor.df[i,1] <- c$estimate
  cor.df[i,3] <- c$p.value
}
for (i in 1:nrow(cor.df)) {
  c <- cor.test(get(paste("w",i,sep=""))[["sum_meanT"]],get(paste("w",i,sep=""))[["GSLthreedayneg3C"]], method = "spearman")
  cor.df[i,2] <- c$estimate
  cor.df[i,4] <- c$p.value
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$T.GDD.P[i] < 0.05) {
    cor.df$T.GDD.Sig[i] <- "Significant"
  } else {
    cor.df$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df)) {
  if (cor.df$T.GSL.P[i] < 0.05) {
    cor.df$T.GSL.Sig[i] <- "Significant"
  } else {
    cor.df$T.GSL.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
cor.df$windowstart <- rep(1982:2009, time = 1)

# Plot
pdf("figs/Correlations10yr_Spearman.pdf", width = 4, height = 3)
ggplot(data = cor.df) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 1988, linetype = "dotted") +
  geom_vline(xintercept = 1997, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GDD.R, shape = T.GDD.Sig),colour = "red",size = 2) +
  geom_line(aes(x = windowstart, y = T.GDD.R),colour = "red",size = 0.5) +
  geom_point(aes(x = windowstart, y = T.GSL.R, shape = T.GSL.Sig),colour = "blue",size = 2) +
  geom_line(aes(x = windowstart, y = T.GSL.R),colour = "blue",size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1991, 2000, 2009)) +
  scale_colour_manual(values = c("blue","red"),
                      labels = c("Temp, GDD","Temp, GSL")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("p > 0.05","p < 0.05")) +
  theme_bw() +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1),
                     limits = c(-0.3, 1)) +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()

# Breakpoint regression of window start year and correlation coefficient
bp.df.10 <- as.data.frame(cor.df$T.GSL.R)
CE.Normal.Mean(bp.df.10, penalty = "AIC", distyp = 1)
CE.Normal.Mean(bp.df.10, penalty = "BIC", distyp = 1)
CE.Normal.Mean(bp.df.10, penalty = "AIC", distyp = 2)
CE.Normal.Mean(bp.df.10, penalty = "BIC", distyp = 2)
# Says there's one break point at 13 (window start 1994) (same result)
CE.Normal.MeanVar(bp.df.10, penalty = "AIC", distyp = 1)
CE.Normal.MeanVar(bp.df.10, penalty = "BIC", distyp = 1)
CE.Normal.MeanVar(bp.df.10, penalty = "AIC", distyp = 2)
CE.Normal.MeanVar(bp.df.10, penalty = "BIC", distyp = 2)
# Says there is one breakpoint at 12 (1993) or three break points at 9, 15, 22 (window start 1990, 1996, 2003)

# Segmented
m <- lm(T.GSL.R ~ windowstart, data = cor.df)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1994)))
summary(my.seg1)
my.seg2 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994)))
summary(my.seg2)
my.seg3 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994, 1998)))
summary(my.seg3)
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 2 breakpoints - 1988, 1997



### 15 yr windows
cor.df.15 <- as.data.frame(matrix(NA, 23, 6))
colnames(cor.df.15) <- c("T.GDD.R","T.GSL.R","T.GDD.P","T.GSL.P","T.GDD.Sig","T.GSL.Sig")
for (i in 1:nrow(cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["sum_GDD"]], method = "spearman")
  cor.df.15[i,1] <- c$estimate
  cor.df.15[i,3] <- c$p.value
}
for (i in 1:nrow(cor.df.15)) {
  c <- cor.test(get(paste("x",i,sep=""))[["sum_meanT"]],get(paste("x",i,sep=""))[["GSLthreedayneg3C"]], method = "spearman")
  cor.df.15[i,2] <- c$estimate
  cor.df.15[i,4] <- c$p.value
}
for (i in 1:nrow(cor.df.15)) {
  if (cor.df.15$T.GDD.P[i] < 0.05) {
    cor.df.15$T.GDD.Sig[i] <- "Significant"
  } else {
    cor.df.15$T.GDD.Sig[i] <- "Not Significant"
  }
}
for (i in 1:nrow(cor.df.15)) {
  if (cor.df.15$T.GSL.P[i] < 0.05) {
    cor.df.15$T.GSL.Sig[i] <- "Significant"
  } else {
    cor.df.15$T.GSL.Sig[i] <- "Not Significant"
  }
}

# Add the windows for the x axis
cor.df.15$windowstart <- rep(1982:2004, time = 1)

# Plot
pdf("figs/Correlations15yr.pdf", width = 4, height = 3)
ggplot(data = cor.df.15) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 1988, linetype = "dotted") +
  geom_vline(xintercept = 1997, linetype = "dotted") +
  geom_point(aes(x = windowstart, y = T.GDD.R, shape = T.GDD.Sig),colour = "red",size = 2) +
  geom_line(aes(x = windowstart, y = T.GDD.R),colour = "red",size = 0.5) +
  geom_point(aes(x = windowstart, y = T.GSL.R, shape = T.GSL.Sig),colour = "blue",size = 2) +
  geom_line(aes(x = windowstart, y = T.GSL.R),colour = "blue",size = 0.5) +
  labs(y = "Correlation Coefficient",
       x = "Window Start Year") +
  scale_x_continuous(breaks = c(1982, 1989, 1996, 2003)) +
  scale_colour_manual(values = c("blue","red"),
                      labels = c("Temp, GDD","Temp, GSL")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("p > 0.05","p < 0.05")) +
  theme_bw() +
  scale_y_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1),
                     limits = c(-0.3, 1)) +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(size = 12))
dev.off()



# Breakpoint regression of window start year and correlation coefficient
bp.df.15 <- as.data.frame(cor.df.15$T.GSL.R)
CE.Normal.Mean(bp.df.15, penalty = "AIC", distyp = 1)
CE.Normal.Mean(bp.df.15, penalty = "BIC", distyp = 1)
CE.Normal.Mean(bp.df.15, penalty = "AIC", distyp = 2)
CE.Normal.Mean(bp.df.15, penalty = "BIC", distyp = 2)
# Says there's one break point at 12 (window start 1994)
CE.Normal.MeanVar(bp.df.15, penalty = "AIC", distyp = 1)
CE.Normal.MeanVar(bp.df.15, penalty = "BIC", distyp = 1)
CE.Normal.MeanVar(bp.df.15, penalty = "AIC", distyp = 2)
CE.Normal.MeanVar(bp.df.15, penalty = "BIC", distyp = 2)
# Says there's one break point at 11 (window start 1993)

# Segmented (I think this is best)
m <- lm(T.GSL.R ~ windowstart, data = cor.df)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1994)))
summary(my.seg1)
my.seg2 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994)))
summary(my.seg2)
my.seg3 <- segmented(m,
                     seg.Z = ~ windowstart, 
                     psi = list(windowstart = c(1992, 1994, 1998)))
summary(my.seg3)
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 2 breakpoints - 1988, 1997



#### Breakpoint PCA ####
# Calculate 3 PCAs for the two time periods before and after the correlation breakpoint
# Use the 10-yr window breakpoints (1990, 1997)
# If D1 part has been run, need to remove those
climateSummer <- select(climateSummer, -d1wp, -d1temp, -d1neg3GSL)
period1 <- subset(climateSummer, year < 1990)
period2 <- subset(climateSummer, year < 1997 & year >= 1990)
period3 <- subset(climateSummer, year >= 1997)

PCAp1 <- rda(na.exclude(period1[,2:ncol(period1)]), scale=T)
p1 <- summary(PCAp1)
p1 # 65.31%, 17.10%
par(mfrow = c(1,1))
biplot(PCAp1, main = "1982-1989")
set.seed(1009)
hornpa(k = 9, size = nrow(period1), reps = 1000, seed = 1009)
eigenvals(PCAp1) # Keep PC1

PCAp2 <- rda(na.exclude(period2[,2:ncol(period2)]), scale=T)
p2 <- summary(PCAp2)
p2 # 48.68%, 24.95%
par(mfrow = c(1,1))
biplot(PCAp2, main = "1990-1996")
set.seed(1009)
hornpa(k = 9, size = nrow(period2), reps = 1000, seed = 1009)
eigenvals(PCAp2) # Keep PC1

PCAp3 <- rda(na.exclude(period3[,2:ncol(period3)]), scale=T)
p3 <- summary(PCAp3)
p3 # 53.66%, 15.08%
par(mfrow = c(1,1))
biplot(PCAp3, main = "1997-2018")
set.seed(1009)
hornpa(k = 9, size = nrow(period3), reps = 1000, seed = 1009)
eigenvals(PCAp3) # Keep PC1

pdf("figs/BreakpointPCAs.pdf", width = 10, height = 8)
par(mfrow = c(1,3))
biplot(PCAp1, xlab = "PC1 65.3%", ylab = "PC2 17.1%")
title("1982-1989", line = 0.2)
biplot(PCAp2, xlab = "PC1 48.7%", ylab = "PC2 = 25.0%")
title("1990 - 1996", line = 0.2)
biplot(PCAp3, xlab = "PC1 53.7%", ylab = "PC2 15.1%")
title("1997 - 2018", line = 0.2)
dev.off()

BP_PCscores <- as.data.frame(rbind(p1$sites[,1:2], p2$sites[,1:2], p3$sites[,1:2])) %>%
  mutate(Year = seq(1982, 2018, 1))

# See how these scores compare with the all years PCA
plot(abs(sumallyrsOutput$sumallPC1), abs(BP_PCscores$PC1))
plot(abs(sumallyrsOutput$sumallPC2), abs(BP_PCscores$PC2))

# So this was an interesting exercise but I am not going to add it to the paper. The PCAs are too different so it doesn't make sense using PC1 from each of them. It is also problematic that two periods have fewer samples than number of variables. This would basically be another way to show the non-stationarity but we have already done that with the moving window PCA and moving window correlations as well as extracted a mean PC score for each year from all windows containing that year. 



############################### Individual Breakpoints #########################
# Breakpoint analysis for each individual climate variable
# Plot each variable
climateSummer_long <- melt(climateSummer,
                           measure.vars = names(climateSummer)[2:10],
                           id.vars = "year")
ggplot(climateSummer_long, aes(year, value)) +
  geom_point() +
  geom_smooth() + 
  labs(x = "Year",
       y = NULL) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y")

# Pettit
library(BreakPoints)
pettit(climateSummer$sum_precip, n_period = 5) # 14, not sig
pettit(climateSummer$sum_moisturedeficit, n_period = 5) # 21, not sig
pettit(climateSummer$fivedayrunning5C, n_period = 5) # 20, not sig
pettit(climateSummer$iceoff_GL4, n_period = 5) # 19, not sig

pettit(climateSummer$sum_meanT, n_period = 5) # 19, sig *
pettit(climateSummer$sum_PET, n_period = 5) # 19, sig *
pettit(climateSummer$sum_GDD, n_period = 5) # 19, sig *
pettit(climateSummer$fivedayrunning12C, n_period = 5) # 19, sig *
pettit(climateSummer$GSLthreedayneg3C, n_period = 5) # 16, sig *

# Segmented
set.seed(923)
m <- lm(sum_meanT ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 2002
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # variable
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # 1988, 1992, 2001
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 3 breakpoints - 1988, 1992, 2001

m <- lm(sum_precip ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 2011
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1987, 1998
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # 1986, 1995, 1997
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(sum_moisturedeficit ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 2014
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1984, 2014
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # None estimated
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(sum_PET ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 2012
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # None estimated
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # None estimated
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(sum_GDD ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 2007
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1986, 2007
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # None estimated
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(fivedayrunning5C ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 1985
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # None estimated
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # None estimated
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(fivedayrunning12C ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 1985.6 or none
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1987, 1992
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # 1987, 1998, 2000 or none
anova(m, my.seg1, my.seg2, my.seg3)
# Best model is 2 breakpoints - 1987, 1992

m <- lm(GSLthreedayneg3C ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 1990
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1993, 1996
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # 1993, 1996, 2012
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

m <- lm(iceoff_GL4 ~ year, data = climateSummer)
summary(m)
my.seg1 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(2000)))
summary(my.seg1) # 1985
my.seg2 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000)))
summary(my.seg2) # 1987, 1993
my.seg3 <- segmented(m,
                     seg.Z = ~ year, 
                     psi = list(year = c(1990, 2000, 2010)))
summary(my.seg3) # 1987, 1997, 2000
anova(m, my.seg1, my.seg2, my.seg3)
# N.S.

# Graph
facet_names <-  c(`sum_meanT` = "Temp",`sum_GDD` = "GDD",`sum_PET` = "PET",`iceoff_GL4` = "Iceoff",`fivedayrunning12C` = "Days 12C",`fivedayrunning5C` = "Days 5C",`sum_precip` = "Precip",`sum_moisturedeficit` = "Moisture Deficit",`GSLthreedayneg3C` = "GSL")
pdf("figs/BreakpointVariables.pdf", width = 8, height = 8)
ggplot(climateSummer_long, aes(year, value)) +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 1988), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 1992), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 2001), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "fivedayrunning12C"),
             aes(xintercept = 1987), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "fivedayrunning12C"),
             aes(xintercept = 1992), linetype = "dotted") +
  geom_point() +
  geom_smooth() + 
  labs(x = "Year",
       y = NULL) +
  facet_wrap(~ variable, ncol = 3, scales = "free", 
             labeller = as_labeller(facet_names))
dev.off()

# Linear graph with breakpoint lines
# Show only significant lines
summary(lm(sum_meanT ~ year, data = climateSummer)) # Sig (break)
summary(lm(sum_precip ~ year, data = climateSummer)) # N.S.
summary(lm(sum_moisturedeficit ~ year, data = climateSummer)) # N.S.
summary(lm(sum_PET ~ year, data = climateSummer)) # Sig
summary(lm(sum_GDD ~ year, data = climateSummer)) # Sig
summary(lm(fivedayrunning5C ~ year, data = climateSummer)) # N.S.
summary(lm(fivedayrunning12C ~ year, data = climateSummer)) # Sig (break)
summary(lm(GSLthreedayneg3C ~ year, data = climateSummer)) # Sig
summary(lm(iceoff_GL4 ~ year, data = climateSummer)) # N.S.

facet_names <-  c(`sum_meanT` = "Temp",`sum_GDD` = "GDD",`sum_PET` = "PET",`iceoff_GL4` = "Iceoff",`fivedayrunning12C` = "Days 12C",`fivedayrunning5C` = "Days 5C",`sum_precip` = "Precip",`sum_moisturedeficit` = "Moisture Deficit",`GSLthreedayneg3C` = "GSL")
pdf("figs/FigureS5.pdf", width = 8, height = 8)
ggplot(climateSummer_long, aes(year, value)) +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 1988), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 1992), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "sum_meanT"),
             aes(xintercept = 2001), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "fivedayrunning12C"),
             aes(xintercept = 1987), linetype = "dotted") +
  geom_vline(data = subset(climateSummer_long, variable == "fivedayrunning12C"),
             aes(xintercept = 1992), linetype = "dotted") +
  geom_point() +
  geom_smooth(data = subset(climateSummer_long, variable == "sum_PET" | 
                              variable == "sum_GDD" | 
                              variable == "GSLthreedayneg3C"),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "sum_meanT" & 
                              year <= 1988),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "sum_meanT" & 
                              year >= 1988 &
                              year <= 1992),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "sum_meanT" & 
                              year >= 1992 &
                              year <= 2001),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "sum_meanT" & 
                              year >= 2001),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "fivedayrunning12C" & 
                              year <= 1987),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "fivedayrunning12C" & 
                              year >= 1987 &
                              year <= 1992),
              method = lm, se = F) + 
  geom_smooth(data = subset(climateSummer_long, variable == "fivedayrunning12C" & 
                              year >= 1992),
              method = lm, se = F) + 
  labs(x = "Year",
       y = NULL) +
  facet_wrap(~ variable, ncol = 3, scales = "free", 
             labeller = as_labeller(facet_names))
dev.off()



#### End Script ####