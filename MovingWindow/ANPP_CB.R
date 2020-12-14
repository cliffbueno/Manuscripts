# Saddle ANPP by PC1 
# 1992 - 2018
# by Cliff Bueno de Mesquita, 2019-2020
# Analysis for NWT advisory visit and non-stationarity manuscript

################################## SETUP #######################################
# clean up enviro, read in needed libraries
rm(list=ls())
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/cliffbueno/Functions/master/Summary.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
library(readxl)
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(r2glmm)
library(reshape2)
library(scales)
library(FSA)
library(corrplot)
library(Rmisc)
library(rsq)
formatter <- function(...){
  function(x) format(round(x, 0), ...)
}
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

#### _FUNCTIONS ####
# set up functions to read in tabular datasets from EDI dynamically
# SCE + CTW code to determine most recent version of package ID and read in current dataset on EDI
# function to determine current version of data package on EDI
getCurrentVersion <- function(edi_id){
   versions=readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-nwt/', edi_id), warn=FALSE)
   currentV <- max(as.numeric(versions))
   return(currentV)
 }
# function to get entity ID for current version
getEntityId <- function(edi_id, version){
   entID <- readLines(paste0('https://pasta.lternet.edu/package/eml/knb-lter-nwt/', edi_id, "/", version, "/"), warn=FALSE)[1]
   entID <- gsub(paste0("http.*/",edi_id,"/",version,"/"), "", entID) # remove all chars except what comes after last /
   return(entID)
 }
# reads in tabular dataset for data package that has only one csv data file (could make more generic with read table, but should know what you're reading in to use)
getTabular <- function(edi_id, na_vals = c("", "NA", NA, NaN, ".", "NaN", " ")){
   v <- getCurrentVersion(edi_id)
   id <- getEntityId(edi_id, v)
   dat <- read.csv(paste0("https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-nwt.", edi_id, ".", v, 
                          "&entityid=", id),
                   strip.white =TRUE, na.strings = na_vals)
   print(paste0("Reading in knb-lter-nwt.", edi_id, ".", v))
   return(dat)
 }
 
source("/Users/cliffbuenodemesquita/Documents/GitHub/long-term-trends/utility_functions/utility_functions_all.R")
source("/Users/cliffbuenodemesquita/Desktop/Functions/r2_nakagawa.R")

####  _DATA ####
# current extended summer PC scores (sumallPC1 = ES PC score)
# For Cliff
extsum <- read.csv("/Users/cliffbuenodemesquita/Documents/GitHub/long-term-trends/extended_summer/NWT_sumallPCclimate_19822018.csv", strip.white = T, na.strings = na_vals)
# For Anna
# extsum <- read.csv("/Users/annawright/Documents/long-term-trends/extended_summer/NWT_sumallPCclimate_19822018.csv", strip.white = T, na.strings = na_vals)

#saddle grid ANPP
sdlprod <- getTabular(16)
sdlprodold <- getTabular(158)
sdlprodold <- subset(sdlprodold, location == "saddle")

# Look at extsum variable correlations
env <- extsum[,c(1:26, 33,38)]
M <- cor(env)
corrplot(M, method = "number", type = "lower", number.cex = 0.3)



####################### ANPP by community type with PC1 ########################
anpp_type_year <- summarySE(sdlprod, measurevar ="NPP",groupvars=c("year","veg_class"))
anpp_type_year <- sdlprod %>%
  group_by(year, veg_class) %>%
  summarize(NPP = mean(NPP), N = n(), se = se(NPP))

# Just look at MM, WM, SB, FF, DM
anpp_type_year_commtype <- subset(anpp_type_year, veg_class == "MM" | veg_class == "WM" | veg_class == "SB" | veg_class == "FF" | veg_class == "DM")

PC1 <- select(extsum, eco_year, sumallPC1, sum_precip, lag_PC1)
colnames(PC1) <- c("Year", "PC1","sum_precip","lag_PC1")

# Add PC1 data to dataframe 
anpp_pc <- left_join(x = PC1, y = anpp_type_year_commtype, by = c("Year" = "year"))

# Basic Linear Regressions
sdlprod_pc <- left_join(x = PC1, y = sdlprod, by = c("Year" = "year"))
dm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "WM"))
summary(wm1) # R2 = 0.03, p = 0.04
sb1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "FF"))
summary(ff1) # R2 = 0.05, p = 0.0058

# Linear Mixed Effects Models (for repeated measure, set plot as random)
dm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_pc, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_pc, veg_class == "MM"))
Anova(mm1) # p = 0.04
r2beta(mm1, partial = FALSE) # R2 = 0.007
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_pc, veg_class == "WM"))
Anova(wm1) # p = 0.01
r2beta(wm1, partial = FALSE) # R2 = 0.04
sb1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_pc, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_pc, veg_class == "FF"))
Anova(ff1) # p = 0.003
r2beta(ff1, partial = FALSE) # R2 = 0.05

# Note big changes with PC1 update! Now, no DM decline, and WM and FF increases!

# Graph
ggplot(anpp_pc, aes(PC1, NPP, colour = veg_class)) +
  geom_point(size = 4, alpha = 0.5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin=NPP-se,ymax=NPP+se),width=0,position = position_dodge(0.1)) +
  geom_smooth(method = lm, se = F, 
              data = subset(anpp_pc, veg_class == "WM" | veg_class == "FF")) +
  labs(x = "Extended Summer (PC1 Axis Score)",
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')')),
       colour = NULL) +
  theme(legend.position = c(0.1, 0.8),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.5, "line"),
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))+
  scale_color_viridis_d()



########################### With Corrected DM and FF ###########################
# Emily Farrer and Lauren Hallett corrected for DM and FF biomass for 1992-1997
# Would need to ask them exactly what they did, but the data is available from the 2016 renewal analyses
# Load their dataset (I put it in the GitHub folder)
setwd("~/Documents/GitHub/long-term-trends/saddle_sppcomp_anpp")
d <- read.csv("NWT_SnowXProdCorrected_EF_LH_2015.csv")

# Need to get the same columns
df1 <- select(d, c(year, plot, class_3, anpp))
df2 <- select(sdlprod, c(year,grid_pt,veg_class,NPP))
# Rename the columns to be the same
colnames(df1) <- c("year","grid_pt","veg_class","NPP")
# Only use non-EDI data for 1992-1997
df1 <- subset(df1, year < 1998)
# Only us EDI data for after 1997
df2 <- subset(df2, year > 1997)
# Now we have two data frames with the same columns and no overlapping years so stack them
df3 <- rbind(df1,df2)
# Now do the same analyses as before

# Summarize by year and veg_class
anpp_type_year <- summarySE(df3, measurevar ="NPP",groupvars=c("year","veg_class"))
anpp_type_year <- df3 %>%
  group_by(year, veg_class) %>%
  summarize(NPP = mean(NPP), N = n(), se = se(NPP))

# View number by veg class
table(anpp_type_year$veg_class)

# Just look at MM, WM, SB, FF, DM
anpp_type_year_commtype <- subset(anpp_type_year, veg_class == "MM" | veg_class == "WM" | veg_class == "SB" | veg_class == "FF" | veg_class == "DM")

PC1 <- select(extsum, eco_year, sumallPC1, sum_precip, lag_PC1)
colnames(PC1) <- c("Year", "PC1","sum_precip","lag_PC1")

# Add PC1 data to dataframe 
anpp_pc <- inner_join(x = PC1, y = anpp_type_year_commtype, by = c("Year" = "year"))

# Regressions on all data
sdlprod_pc <- left_join(x = PC1, y = df3, by = c("Year" = "year"))
dm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "WM"))
summary(wm1) # R2 = 0.03, p = 0.04
sb1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ PC1, data = subset(sdlprod_pc, veg_class == "FF"))
summary(ff1) # NS

# Test adding in summer precip. With interaction.
# see how PC1 and summer precip relate
plot(sdlprod_pc$PC1, sdlprod_pc$sum_precip)
# Not too correlated
dm2 <- lm(NPP ~ PC1*sum_precip, data = subset(sdlprod_pc, veg_class == "DM"))
summary(dm2) # PC1***, Precip***, Interaction***
mm2 <- lm(NPP ~ PC1*sum_precip, data = subset(sdlprod_pc, veg_class == "MM"))
summary(mm2) # Precip**
wm2 <- lm(NPP ~ PC1*sum_precip, data = subset(sdlprod_pc, veg_class == "WM"))
summary(wm2) # PC1*, Interaction*
sb2 <- lm(NPP ~ PC1*sum_precip, data = subset(sdlprod_pc, veg_class == "SB"))
summary(sb2) # Precip***
ff2 <- lm(NPP ~ PC1*sum_precip, data = subset(sdlprod_pc, veg_class == "FF"))
summary(ff2) # NS

# Now no interaction.
dm3 <- lm(NPP ~ PC1+sum_precip, data = subset(sdlprod_pc, veg_class == "DM"))
summary(dm3) # Precip*** (negative)
mm3 <- lm(NPP ~ PC1+sum_precip, data = subset(sdlprod_pc, veg_class == "MM"))
summary(mm3) # Precip*** (negative)
wm3 <- lm(NPP ~ PC1+sum_precip, data = subset(sdlprod_pc, veg_class == "WM"))
summary(wm3) # PC1.(marginally negative)
sb3 <- lm(NPP ~ PC1+sum_precip, data = subset(sdlprod_pc, veg_class == "SB"))
summary(sb3) # Precip** (negative)
ff3 <- lm(NPP ~ PC1+sum_precip, data = subset(sdlprod_pc, veg_class == "FF"))
summary(ff3) # Precip* (negative)

# Test for lag
dm4 <- lm(NPP ~ lag_PC1, data = subset(sdlprod_pc, veg_class == "DM"))
summary(dm4) # R2 = 0.02, p = 0.008
mm4 <- lm(NPP ~ lag_PC1, data = subset(sdlprod_pc, veg_class == "MM"))
summary(mm4) # NS
wm4 <- lm(NPP ~ lag_PC1, data = subset(sdlprod_pc, veg_class == "WM"))
summary(wm4) # NS
sb4 <- lm(NPP ~ lag_PC1, data = subset(sdlprod_pc, veg_class == "SB"))
summary(sb4) # NS
ff4 <- lm(NPP ~ lag_PC1, data = subset(sdlprod_pc, veg_class == "FF"))
summary(ff4) # NS

# Graph showing increase in WM anpp
ggplot(anpp_pc, aes(PC1, NPP, colour = veg_class)) +
  geom_point(size = 4, alpha = 0.5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin=NPP-se,ymax=NPP+se),width=0,position = position_dodge(0.1)) +
  geom_smooth(method = lm, se = F, data = subset(anpp_pc, veg_class == "WM")) +
  labs(x = "Extended Summer (PC1 Axis Score)",
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')')),
       colour = NULL) +
  theme(legend.position = c(0.1, 0.8),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.25, "line"),
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))

# Graph with summer precipition
ggplot(anpp_pc, aes(sum_precip, NPP, colour = veg_class)) +
  geom_point(size = 4, alpha = 0.5, position = position_dodge(5)) +
  geom_errorbar(aes(ymin=NPP-se,ymax=NPP+se),width=0,position = position_dodge(5)) +
  geom_smooth(method = lm, se = F) +
  labs(x = "Summer Precipitation",
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')')),
       colour = NULL) +
  theme(legend.position = c(0.85, 0.8),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.25, "line"),
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))+
  scale_color_viridis_d()

# Lag Graph showing increase in DM anpp
ggplot(anpp_pc, aes(lag_PC1, NPP, colour = veg_class)) +
  geom_point(size = 4, alpha = 0.5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin=NPP-se,ymax=NPP+se),width=0,position = position_dodge(0.1)) +
  geom_smooth(method = lm, se = F, data = subset(anpp_pc, veg_class == "DM")) +
  labs(x = "Previous Year PC1 Axis Score",
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')')),
       colour = NULL) +
  theme(legend.position = c(0.1, 0.82),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.25, "line"),
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))

# Graph to show interactions between PC1 and Precip for DM, WM
dmwm <- subset(anpp_pc, veg_class == "DM" | veg_class == "WM")
ggplot(dmwm, aes(PC1, NPP, shape = veg_class)) +
  geom_point(aes(col = sum_precip),size = 4, alpha = 0.5, position = position_dodge(0.1)) +
  geom_errorbar(aes(ymin=NPP-se,ymax=NPP+se),width=0,position = position_dodge(0.1)) +
  geom_smooth(se = F) +
  labs(x = "Extended Summer (PC1 Axis Score)",
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  scale_color_viridis_c(name = "Summer Precip", option = "B", direction = -1) +
  theme(legend.position = "right",
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.key.size = unit(1.25, "line"),
        axis.text = element_text(size = 16), 
        axis.title = element_text(face="bold",size=18))



######################### CB Analysis for PCA methods paper ####################
# Test PC1, Temp, GDD, GSL as predictors of biomass. Then test certain interactions. Graph.
# All subsets
# Add vif, gls, randomforests

#### _Setup ####
d <- read.csv("NWT_SnowXProdCorrected_EF_LH_2015.csv")

# Need to get the same columns
df1 <- select(d, c(year, plot, class_3, anpp))
df2 <- select(sdlprod, c(year,grid_pt,veg_class,NPP))
# Rename the columns to be the same
colnames(df1) <- c("year","grid_pt","veg_class","NPP")
# Only use non-EDI data for 1992-1997
df1 <- subset(df1, year < 1998)
# Only us EDI data for after 1997
df2 <- subset(df2, year > 1997)
# Now we have two data frames with the same columns and no overlapping years so stack them
df3 <- rbind(df1,df2)

# Summarize by year and veg_class
anpp_type_year <- summarySE(df3, 
                            measurevar ="NPP",
                            groupvars=c("year","veg_class"))
anpp_type_year <- df3 %>%
  group_by(year, veg_class) %>%
  summarize(NPP = mean(NPP), N = n(), se = se(NPP))

# View number by veg class
table(anpp_type_year$veg_class)

# Just look at MM, WM, SB, FF, DM
anpp_type_year_commtype <- subset(anpp_type_year, veg_class == "MM" | veg_class == "WM" | veg_class == "SB" | veg_class == "FF" | veg_class == "DM")

clim <- select(extsum, c(eco_year, sumallPC1, sum_meanT, sum_GDD, GSLthreedayneg3C, sum_precip, moisturedeficit))
colnames(clim) <- c("Year", "PC1","Temp","GDD","GSL","Precip","MD")

# Add clim data to commtype summary dataframe 
anpp_clim <- inner_join(x = clim, y = anpp_type_year_commtype, by = c("Year" = "year"))

# Add clim data to all plot dataframe
sdlprod_clim <- left_join(x = clim, y = df3, by = c("Year" = "year"))

# Add all extended summer data to all plot dataframe
sdlprod_clim_tot <- left_join(x = extsum, y = df3, by = c("eco_year" = "year"))

#### _All Subsets Regression for WM ####
sdlprod_clim_tot_WM <- subset(sdlprod_clim_tot, veg_class == "WM")
library(bestglm)
env <- select(sdlprod_clim_tot_WM, c(sumallPC1, sumallPC2, sum_meanT, sum_GDD, fivedayrunning5C, fivedayrunning12C,GSLthreedayneg3C,sum_precip,sum_moisturedeficit,sum_PET,iceoff_GL4))
M <- cor(env)
corrplot(M, method = "number", type = "lower")
# Remove GDD and PET based on VIF
# Remove 1997 as below
sdlprod_clim_tot_WM_no1997 <- subset(sdlprod_clim_tot_WM, eco_year != 1997)
env <- select(sdlprod_clim_tot_WM_no1997, c(sumallPC1, sumallPC2, sum_meanT, fivedayrunning5C, fivedayrunning12C,GSLthreedayneg3C,sum_precip,sum_moisturedeficit,iceoff_GL4))
y <- sdlprod_clim_tot_WM_no1997$NPP
Xy <- as.data.frame(cbind(env,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel # Same (GSL and MD)
bestmodel$BestModels
nullmodel <- lm(NPP ~ 1, data = sdlprod_clim_tot_WM)
bestmodel <- lm(NPP ~ GSLthreedayneg3C + sum_moisturedeficit, 
                data = sdlprod_clim_tot_WM)
summary(bestmodel)
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
# Partial R2
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)



#### _All subsets regression revised ####
# Run for two time periods (the first won't have enough data though)
# Use original variables, not axes scores
# Compare all top models - point to show that suite of variables at play
# Remove GDD and PET to keep VIFs < 7.5
period1 <- subset(sdlprod_clim_tot_WM_no1997, eco_year < 2000)
period2 <- subset(sdlprod_clim_tot_WM_no1997, eco_year > 2000)
period1env <- select(period1, 
                     c(sum_meanT, fivedayrunning5C, fivedayrunning12C,
                       GSLthreedayneg3C, sum_precip, sum_moisturedeficit, 
                       iceoff_GL4))
period2env <- select(period2, 
                     c(sum_meanT, fivedayrunning5C, fivedayrunning12C,
                       GSLthreedayneg3C, sum_precip, sum_moisturedeficit, 
                       iceoff_GL4))

# Period 1
y <- period1$NPP
Xy <- as.data.frame(cbind(period1env, y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, 
                     TopModels = 10)
bestmodel # GSL
bestmodel$BestModels
p1bestmods <- as.data.frame(bestmodel$BestModels) %>%
  mutate(Period = "1992 - 1995")
nullmodel <- lm(NPP ~ 1, data = period1)
bestmodel <- lm(NPP ~ GSLthreedayneg3C, data = period1)
summary(bestmodel)
par(mfrow = c(2,2))
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel) # Slight improvement!
anova(nullmodel, bestmodel) # Not significant! (marginal p = 0.06)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)

# Period 2
y <- period2$NPP
Xy <- as.data.frame(cbind(period2env, y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, 
                     TopModels = 10)
bestmodel # Moisture deficit
bestmodel$BestModels # All similar
p2bestmods <- as.data.frame(bestmodel$BestModels) %>%
  mutate(Period = "2008 - 2018")
nullmodel <- lm(NPP ~ 1, data = period2)
bestmodel <- lm(NPP ~ sum_moisturedeficit, data = period2)
summary(bestmodel)
par(mfrow = c(2,2))
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel) # Improvement!
anova(nullmodel, bestmodel) # Significant!
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)

bestmodsup <- rbind(p1bestmods, p2bestmods) %>%
  mutate(Model = rep(1:10, 2),
         AIC = Criterion) %>%
  select(Period, Model, everything(), -Criterion)
write.csv(bestmodsup, "BestModels.csv")



#### _lmer ####
# PC1
dm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # p = 0.04
r2beta(mm1, partial = FALSE) # R2 = 0.007
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # p = 0.013
r2beta(wm1, partial = FALSE) # R2 = 0.04
sb1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# Temp
dm1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # NS
wm1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # p = 0.04916
r2beta(wm1, partial = FALSE) # R2 = 0.025
sb1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# GDD
dm1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # NS
wm1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # p = 0.1166
r2beta(wm1, partial = FALSE) # R2 = 0.016
sb1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# GSL
dm1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # NS
wm1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # p = 0.012
r2beta(wm1, partial = FALSE) # R2 = 0.04
sb1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# Precip
dm1 <- lmer(NPP ~ Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # Sig
mm1 <- lmer(NPP ~ Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # Sig
wm1 <- lmer(NPP ~ Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # p = 0.05752
r2beta(wm1, partial = FALSE) # R2 = 0.023
sb1 <- lmer(NPP ~ Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # Sig
ff1 <- lmer(NPP ~ Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # Sig

# MD
dm1 <- lmer(NPP ~ MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # Sig
mm1 <- lmer(NPP ~ MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # NS
wm1 <- lmer(NPP ~ MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # NS
sb1 <- lmer(NPP ~ MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # Sig
ff1 <- lmer(NPP ~ MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# Interactions
dm1 <- lmer(NPP ~ Temp*GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ Temp*GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # NS
wm1 <- lmer(NPP ~ Temp*GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # NS
sb1 <- lmer(NPP ~ Temp*GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ Temp*GSL + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

dm1 <- lmer(NPP ~ Temp*Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # Precip and Interaction Significant
mm1 <- lmer(NPP ~ Temp*Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # Precip and Interaction Significant
wm1 <- lmer(NPP ~ Temp*Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # All significant
sb1 <- lmer(NPP ~ Temp*Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # Precip significant
ff1 <- lmer(NPP ~ Temp*Precip + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # Precip significant

dm1 <- lmer(NPP ~ Temp*MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "DM"))
Anova(dm1) # MD and Interaction Significant
mm1 <- lmer(NPP ~ Temp*MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "MM"))
Anova(mm1) # Interaction Significant
wm1 <- lmer(NPP ~ Temp*MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "WM"))
Anova(wm1) # NS
sb1 <- lmer(NPP ~ Temp*MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "SB"))
Anova(sb1) # MD significant
ff1 <- lmer(NPP ~ Temp*MD + (1|grid_pt), data = subset(sdlprod_clim, veg_class == "FF"))
Anova(ff1) # NS

# Summarize Wet Meadow. Use this for paper. PC1 and GSL were significant, Temp and GDD were not.
sdlprod_clim_wm <- subset(sdlprod_clim, veg_class == "WM")
sumWM <- ddply(sdlprod_clim_wm, c("Year","PC1","Temp","GDD","GSL"), summarise,
             meanNPP = mean(NPP),
             seNPP = se(NPP))
sdlprod_clim_wm$n <- 1
ncheck <- aggregate(sdlprod_clim_wm$n, by = list(sdlprod_clim_wm$Year), sum)
ncheck

# Facet wrap, 4 columns.
WMfacet <- melt(sumWM, measure.vars = c("PC1","Temp","GDD","GSL"),
                id.vars = c("Year","meanNPP","seNPP"))
facet_names <- c('PC1' = "a) PC1", 'Temp' = "b) Temp", 'GDD' = "c) GDD", 'GSL' = "d) GSL")
ggplot(WMfacet, aes(value, meanNPP)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0, alpha = 0.5) +
  geom_smooth(data = subset(WMfacet, variable == "PC1" | variable == "GSL"), method = lm) +
  facet_wrap(~ variable, ncol = 4, scales = "free_x", labeller = as_labeller(facet_names)) +
  ylab(bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  xlab(NULL) +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size=16),
        strip.text = element_text(size = 14),
        panel.margin = unit(0.95, "lines"))

# 1997 has huge error bars. Only 2 plots, and very different. Remove. Use this!
sumWM2 <- subset(sumWM, Year != 1997)
WMfacet2 <- melt(sumWM2, measure.vars = c("PC1","Temp","GDD","GSL"),
                id.vars = c("Year","meanNPP","seNPP"))
pdf("ANPP_PC_T_GDD_GSL.pdf", width = 7, height = 3)
ggplot(WMfacet2, aes(value, meanNPP)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0, alpha = 0.5) +
  geom_smooth(data = subset(WMfacet2, variable == "PC1" | variable == "GSL"), method = lm) +
  facet_wrap(~ variable, ncol = 4, scales = "free_x", labeller = as_labeller(facet_names)) +
  ylab(bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  xlab(NULL) +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size=16),
        strip.text = element_text(size = 14),
        panel.spacing = unit(0.95, "lines"))
dev.off()

# Test removing 1997 (Results same - use this!)
# Get marginal and conditional R2
library(insight)
library(performance)
sdlprod_clim_no1997 <- subset(sdlprod_clim, Year != 1997)
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim_no1997, veg_class == "WM"))
Anova(wm1) # R2 = 0.04, p = 0.01
r2beta(wm1, partial = FALSE)
r2_nakagawa(wm1) # 0.31, 0.04
wm1 <- lmer(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim_no1997, veg_class == "WM"))
Anova(wm1) # R2 = 0.03, p = 0.05
r2beta(wm1, partial = FALSE)
r2_nakagawa(wm1) # 0.30, 0.02
wm1 <- lmer(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim_no1997, veg_class == "WM"))
Anova(wm1) # R2 = 0.02, p = 0.12
r2beta(wm1, partial = FALSE)
r2_nakagawa(wm1) # 0.29, 0.01
wm1 <- lmer(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim_no1997, veg_class == "WM"))
Anova(wm1) # R2 = 0.04, p = 0.01
r2beta(wm1, partial = FALSE)
r2_nakagawa(wm1) # 0.31, 0.04



#### _Idea 1: Moving window regressions with biomass ####
# Actually can't really do because only 15 years of biomass data (14 if removing 1997)
# So, can just do first half and second half
sdlprod_clim_no1997_wm <- subset(sdlprod_clim_no1997, veg_class == "WM")
oldyears <- subset(sdlprod_clim_no1997_wm, Year < 2000)
newyears <- subset(sdlprod_clim_no1997_wm, Year > 2000)
wmO <- lmer(NPP ~ PC1 + (1|grid_pt), data = oldyears)
Anova(wmO) # p = 0.15
r2beta(wmO, partial = FALSE) # R2 = 0.07
r2_nakagawa(wmO)
wmO <- lmer(NPP ~ Temp + (1|grid_pt), data = oldyears)
Anova(wmO) # p = 0.02*
r2beta(wmO, partial = FALSE) # R2 = 0.16
r2_nakagawa(wmO)
wmO <- lmer(NPP ~ GDD + (1|grid_pt), data = oldyears)
Anova(wmO) # p = 0.24
r2beta(wmO, partial = FALSE) # R2 = 0.05
wmO <- lmer(NPP ~ GSL + (1|grid_pt), data = oldyears)
Anova(wmO) # p = 0.01*
r2beta(wmO, partial = FALSE) # R2 = 0.17
r2_nakagawa(wmO)
# Temp and GSL significant, PC1 and GDD not.

wmO <- lmer(NPP ~ PC1 + (1|grid_pt), data = newyears)
Anova(wmO) # p = 0.10
r2beta(wmO, partial = FALSE) # R2 = 0.02
wmO <- lmer(NPP ~ Temp + (1|grid_pt), data = newyears)
Anova(wmO) # p = 0.74
r2beta(wmO, partial = FALSE) # R2 = 0.001
wmO <- lmer(NPP ~ GDD + (1|grid_pt), data = newyears)
Anova(wmO) # p = 0.87
r2beta(wmO, partial = FALSE) # R2 = 0
wmO <- lmer(NPP ~ GSL + (1|grid_pt), data = newyears)
Anova(wmO) # p = 0.38
r2beta(wmO, partial = FALSE) # R2 = 0.006
# Nothing significant. This is interesting. The old years are driving the result. This is before temp and GSL were decoupled.
# Plot showing old vs new.
WMfacet2$oldnew <- NA
for (i in 1:nrow(WMfacet2)) {
  ifelse(WMfacet2$Year[i] < 2000,
         WMfacet2$oldnew[i] <- "1992-1995",
         WMfacet2$oldnew[i] <- "2008-2018")
}

library(ggrepel)
pdf("figs/ANPP_PC_T_GDD_GSL_OldNew.pdf", width = 7, height = 3)
ggplot(WMfacet2, aes(value, meanNPP, colour = oldnew)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0, alpha = 0.5) +
  geom_smooth(data = subset(WMfacet2, oldnew == "1992-1995" & variable == "Temp"), 
              method = lm, se = FALSE, show.legend = FALSE, size = 0.25) +
  geom_smooth(data = subset(WMfacet2, oldnew == "1992-1995" & variable == "GSL"),
              method = lm, se = FALSE, show.legend = FALSE, size = 0.25) +
  geom_text(data = subset(WMfacet2, variable != "PC1"),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 1992),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 1993),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 1994),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 1995),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 2008),
            aes(label = Year), size = 2, color = "black", nudge_x = 0.06, nudge_y = 1.25) +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 2010),
            aes(label = Year), size = 2, color = "black") +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year == 2011),
            aes(label = Year), size = 2, color = "black", nudge_x = -0.06, nudge_y = -1.25) +
  geom_text(data = subset(WMfacet2, variable == "PC1" & Year > 2011),
            aes(label = Year), size = 2, color = "black") +
  scale_x_continuous(expand = c(0.08, 0.08)) +
  facet_wrap(~ variable, ncol = 4, scales = "free_x", labeller = as_labeller(facet_names)) +
  ylab(bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  xlab(NULL) +
  labs(colour = NULL) +
  scale_colour_manual(values = c("blue","red")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position = c(0.92, 0.1),
        legend.background = element_rect(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.05, "cm"),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size=16),
        strip.text = element_text(size = 14),
        panel.spacing = unit(0.5, "lines"))
dev.off()

# Update with color for each year, shape by the two periods
# Looks terrible don't use.
pdf("Figure5.pdf", width = 7, height = 3)
ggplot(WMfacet2, aes(value, meanNPP, colour = Year, shape = oldnew)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0, alpha = 0.5) +
  geom_smooth(data = subset(WMfacet2, oldnew == "1992-1995" & variable == "Temp"), 
              method = lm, se = FALSE, show.legend = FALSE) +
  geom_smooth(data = subset(WMfacet2, oldnew == "1992-1995" & variable == "GSL"),
              method = lm, se = FALSE, show.legend = FALSE) +
  facet_wrap(~ variable, ncol = 4, scales = "free_x", labeller = as_labeller(facet_names)) +
  ylab(bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  xlab(NULL) +
  labs(colour = NULL) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2000) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  theme(legend.position = "right",
        legend.background = element_rect(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.05, "cm"),
        axis.text = element_text(size = 10), 
        axis.title = element_text(size=16),
        strip.text = element_text(size = 14),
        panel.spacing = unit(0.95, "lines"))
dev.off()



#### _Idea 2: Test Mean Moving Window PC Scores ####
# This was an idea for how to incorporate non-stationarity into the analysis
# Data (produced by MovingWindow.R script)
mean_mw_scores <- read.csv("../extended_summer/sensitivity_analyses/mean_window_scores.csv",
                           row.names = 1)
# Merge with biomass
df <- select(sdlprod_clim, -PC1)
df <- left_join(df, mean_mw_scores, by = "Year")

# Test new PC1
# PC1
dm1 <- lm(NPP ~ PC1, data = subset(df, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ PC1, data = subset(df, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ PC1, data = subset(df, veg_class == "WM"))
summary(wm1) # R2 = 0.03, p = 0.078
sb1 <- lm(NPP ~ PC1, data = subset(df, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ PC1, data = subset(df, veg_class == "FF"))
summary(ff1) # R2 = 0.02, p = 0.0361

dm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(df, veg_class == "DM"))
Anova(dm1) # NS
mm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(df, veg_class == "MM"))
Anova(mm1) # p = 0.001522
r2beta(mm1, partial = FALSE) # R2 = 0.017
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(df, veg_class == "WM"))
Anova(wm1) # p = 0.03351
r2beta(wm1, partial = FALSE) # R2 = 0.029
sb1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(df, veg_class == "SB"))
Anova(sb1) # NS
ff1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(df, veg_class == "FF"))
Anova(ff1) # p = 0.02614
r2beta(ff1, partial = FALSE) # R2 = 0.23

# Check WM no 1997
dfno97 <- subset(df, Year != 1997)
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = subset(dfno97, veg_class == "WM"))
Anova(wm1) # p = 0.03323
r2beta(wm1, partial = FALSE) # R2 = 0.03



#### _Revisions: vif, gls, random forests ####
# VIF 
m_allvars <- lm(sdlprod_clim_tot$NPP ~ sdlprod_clim_tot$sum_meanT + sdlprod_clim_tot$sum_GDD + sdlprod_clim_tot$fivedayrunning5C + sdlprod_clim_tot$fivedayrunning12C + sdlprod_clim_tot$GSLthreedayneg3C + sdlprod_clim_tot$sum_precip + sdlprod_clim_tot$sum_moisturedeficit + sdlprod_clim_tot$sum_PET + sdlprod_clim_tot$iceoff_GL4)
vif(m_allvars)
m_sub1 <- lm(sdlprod_clim_tot$NPP ~ sdlprod_clim_tot$sum_meanT + sdlprod_clim_tot$fivedayrunning5C + sdlprod_clim_tot$fivedayrunning12C + sdlprod_clim_tot$GSLthreedayneg3C + sdlprod_clim_tot$sum_precip + sdlprod_clim_tot$sum_moisturedeficit + sdlprod_clim_tot$sum_PET + sdlprod_clim_tot$iceoff_GL4)
vif(m_sub1)
m_sub2 <- lm(sdlprod_clim_tot$NPP ~ sdlprod_clim_tot$sum_meanT + sdlprod_clim_tot$fivedayrunning5C + sdlprod_clim_tot$fivedayrunning12C + sdlprod_clim_tot$GSLthreedayneg3C + sdlprod_clim_tot$sum_precip + sdlprod_clim_tot$sum_moisturedeficit + sdlprod_clim_tot$iceoff_GL4)
vif(m_sub2) # Removing GDD and PET keeps VIFs below 7.5

# GLS correlation
# - GLS doesn't accept random factors, so trying correlation argument in lme.
sdlprod_clim_no1997_wm <- subset(sdlprod_clim_no1997, veg_class == "WM")
library(nlme)
gls1 <- lme(NPP ~ PC1,
            random = ~1|grid_pt, 
            data = sdlprod_clim_no1997_wm,
            method = "REML")
Anova(gls1) # R2 = 0.04, p = 0.01 (this is same as lmer result)
gls1.cor <- lme(NPP ~ PC1,
            random = ~1|grid_pt,
            correlation = corAR1(form = ~ Year),
            data = sdlprod_clim_no1997_wm,
            method = "REML") # Error

# Random Forests
library(randomForest)
rf <- randomForest(NPP ~ sumallPC1 + sumallPC2 + sum_meanT + fivedayrunning5C + fivedayrunning12C + GSLthreedayneg3C + sum_precip + sum_moisturedeficit + iceoff_GL4, 
                   data = sdlprod_clim_tot_WM_no1997, 
                   ntree = 5000, 
                   importance = TRUE)
rf
importance(rf)
varImpPlot(rf)
# Try with random effect
library(MixRF)
rf.mm <- MixRF(Y = sdlprod_clim_tot_WM_no1997$NPP, 
               X = env,
               random = "(1|grid_pt)", 
               data = sdlprod_clim_tot_WM_no1997, 
               initialRandomEffects = 0, 
               ErrorTolerance = 0.001, 
               MaxIterations = 1000)
rf.mm
summary(rf.mm)
importance(rf.mm)
varImpPlot(rf)

### _tests/explore ####
# Try getting the sum ANPP for each year.
summedNPP <- ddply(sdlprod_clim, c("Temp","Precip"), summarise,
                   summedNPP = sum(NPP))
m1 <- lm(summedNPP ~ Temp*Precip, data = summedNPP)
summary(m1) # NS

# Summarize Dry Meadow
library(FSA)
sdlprod_clim_dm <- subset(sdlprod_clim, veg_class == "DM")
sum <- ddply(sdlprod_clim_dm, c("Temp","Precip"), summarise,
             meanNPP = mean(NPP),
             seNPP = se(NPP))

# Temperature
ggplot(sum, aes(Temp, meanNPP)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0) +
  labs(x = expression(bold("Mean Summer Temperature " ( degree*C))),
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18))

# Precip
ggplot(sum, aes(Precip, meanNPP)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0) +
  labs(x = expression(bold("Summer Precipitation (mm)")),
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18))

# Interaction
ggplot(sum, aes(Temp, meanNPP, colour = Precip)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin=meanNPP-seNPP,ymax=meanNPP+seNPP),width=0) +
  geom_smooth(data = subset(sum, Precip < 150), method = lm, se = F, color = "red") +
  geom_smooth(data = subset(sum, Precip > 250), method = lm, se = F, color = "blue") +
  geom_smooth(data = subset(sum, Precip < 250&Precip > 150), method = lm, se = F, color = "green") +
  labs(x = expression(bold("Mean Summer Temperature " ( degree*C))),
       y = bquote(bold('Mean ANPP ('*g~m^-2~yr^-1*')'))) +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size=18))

# lms
dm1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # R2 = 0.03, p = 0.04
sb1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # R2 - 0.04, P = 0.03
sb1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ Precip, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # R2 = 0.08, p < 0.001
mm1 <- lm(NPP ~ Precip, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # R2 = 0.02, p < 0.001
wm1 <- lm(NPP ~ Precip, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Precip, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # R2 = 0.02, p < 0.01
ff1 <- lm(NPP ~ Precip, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # R2 = 0.02, p = 0.03
dm1 <- lm(NPP ~ MD, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # R2 = 0.01, p = 0.04
mm1 <- lm(NPP ~ MD, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ MD, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ MD, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # R2 = 0.01, p = 0.02
ff1 <- lm(NPP ~ MD, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ Temp*GSL, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ Temp*GSL, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp*GSL, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp*GSL, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*GSL, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ Temp*Precip, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # Temp, Precip and Interaction Significant
mm1 <- lm(NPP ~ Temp*Precip, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # Marginal interaction
wm1 <- lm(NPP ~ Temp*Precip, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # Temp significant, interaction marginal
sb1 <- lm(NPP ~ Temp*Precip, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*Precip, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS
dm1 <- lm(NPP ~ Temp*MD, data = subset(sdlprod_clim, veg_class == "DM"))
summary(dm1) # Significant interaction
mm1 <- lm(NPP ~ Temp*MD, data = subset(sdlprod_clim, veg_class == "MM"))
summary(mm1) # Marginal interaction
wm1 <- lm(NPP ~ Temp*MD, data = subset(sdlprod_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp*MD, data = subset(sdlprod_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*MD, data = subset(sdlprod_clim, veg_class == "FF"))
summary(ff1) # NS

# Test polynomial. None significant
m1 <- lm(sdlprod_clim_wm$NPP ~ poly(sdlprod_clim_wm$PC1,2,raw=TRUE))
summary(m1)
m1 <- lm(sdlprod_clim_wm$NPP ~ poly(sdlprod_clim_wm$Temp,2,raw=TRUE))
summary(m1)
m1 <- lm(sdlprod_clim_wm$NPP ~ poly(sdlprod_clim_wm$GDD,2,raw=TRUE))
summary(m1)
m1 <- lm(sdlprod_clim_wm$NPP ~ poly(sdlprod_clim_wm$GSL,2,raw=TRUE))
summary(m1)

# Test removing 2012 - decide to keep!
sdlprod_clim_no2012 <- subset(sdlprod_clim, Year != 2012)
wm1 <- lm(NPP ~ PC1, data = subset(sdlprod_clim_no2012, veg_class == "WM"))
summary(wm1) # R2 = 0.04, p = 0.02
wm1 <- lm(NPP ~ Temp, data = subset(sdlprod_clim_no2012, veg_class == "WM"))
summary(wm1) # R2 = 0.03, p = 0.07
wm1 <- lm(NPP ~ GDD, data = subset(sdlprod_clim_no2012, veg_class == "WM"))
summary(wm1) # R2 = 0.02, p = 0.15
wm1 <- lm(NPP ~ GSL, data = subset(sdlprod_clim_no2012, veg_class == "WM"))
summary(wm1) # R2 = 0.04, p = 0.03
# So, the results are the same!
wm1 <- lm(NPP ~ PC1 + (1|grid_pt), data = subset(sdlprod_clim_no2012, veg_class == "WM"))
Anova(wm1) # R2 = 0.04, p = 0.02
r2beta(wm1, partial = FALSE)
wm1 <- lm(NPP ~ Temp + (1|grid_pt), data = subset(sdlprod_clim_no2012, veg_class == "WM"))
Anova(wm1) # R2 = 0.03, p = 0.08
r2beta(wm1, partial = FALSE)
wm1 <- lm(NPP ~ GDD + (1|grid_pt), data = subset(sdlprod_clim_no2012, veg_class == "WM"))
Anova(wm1) # R2 = 0.02, p = 0.15
r2beta(wm1, partial = FALSE)
wm1 <- lm(NPP ~ GSL + (1|grid_pt), data = subset(sdlprod_clim_no2012, veg_class == "WM"))
Anova(wm1) # R2 = 0.04, p = 0.03
r2beta(wm1, partial = FALSE)
# Results same. R2 computation not working though.

# Starting in 2008, there are two clipping plots per plot. Get average.
sdlprod_clim_no1997_wm <- subset(sdlprod_clim_no1997, veg_class == "WM")
sdlprod_clim_no1997_wm_avg <- ddply(sdlprod_clim_no1997_wm, 
                                    c("Year","PC1","Temp","GDD","GSL","grid_pt"), summarise,
                                    NPP = mean(NPP))
wm1 <- lm(NPP ~ PC1, data = sdlprod_clim_no1997_wm_avg)
summary(wm1) # R2 = 0.06, p = 0.0358
wm1 <- lm(NPP ~ Temp, data = sdlprod_clim_no1997_wm_avg)
summary(wm1) # R2 = 0.05, p = 0.06
wm1 <- lm(NPP ~ GDD, data = sdlprod_clim_no1997_wm_avg)
summary(wm1) # R2 = 0.03, p = 0.13
wm1 <- lm(NPP ~ GSL, data = sdlprod_clim_no1997_wm_avg)
summary(wm1) # R2 = 0.08, p = 0.0176
# So, the results are the same!
wm1 <- lmer(NPP ~ PC1 + (1|grid_pt), data = sdlprod_clim_no1997_wm_avg)
Anova(wm1) # R2 = 0.07, p = 0.00879
r2beta(wm1, partial = FALSE)
wm1 <- lmer(NPP ~ Temp + (1|grid_pt), data = sdlprod_clim_no1997_wm_avg)
Anova(wm1) # R2 = 0.05, p = 0.02
r2beta(wm1, partial = FALSE)
wm1 <- lmer(NPP ~ GDD + (1|grid_pt), data = sdlprod_clim_no1997_wm_avg)
Anova(wm1) # R2 = 0.03, p = 0.07
r2beta(wm1, partial = FALSE)
wm1 <- lmer(NPP ~ GSL + (1|grid_pt), data = sdlprod_clim_no1997_wm_avg)
Anova(wm1) # R2 = 0.08, p = 0.002759
r2beta(wm1, partial = FALSE)
# Results basically the same except temp. I think it's better to use both plots, and there is already the random factor of plot.



################################## Other #######################################
# Realized that on EDI the data ID 16 only goes back to 1992
# ID 158 has 1982-2011. So combine those two? But what is ID 158??
# If using WM then don't need to worry about correcting the data
# Need to get the same columns
df1 <- select(sdlprod, c(year,grid_pt,veg_class,NPP))
df2 <- select(sdlprodold, c(year,plotid,veg_class,NPP))
# Rename the columns to be the same
colnames(df2) <- c("year","grid_pt","veg_class","NPP")
# Only use EDI 158 data for 1982-2011
df1 <- subset(df1, year > 2011)
# Only use EDI 16 data for after 2011
df2 <- subset(df2, year < 2012)
# Now we have two data frames with the same columns and no overlapping years so stack them
df3 <- rbind(df1,df2)

# Summarize by year and veg_class and plot id (some years had subsamples)
anpp_type_year <- summarySE(df3, measurevar ="NPP",
                            groupvars=c("year","veg_class","grid_pt"))

# View number by veg class
table(anpp_type_year$veg_class)

# Just look at MM, WM, SB, FF
anpp_type_year_commtype <- subset(anpp_type_year, veg_class == "MM" | veg_class == "WM" | veg_class == "SB" | veg_class == "FF")

clim <- select(extsum, c(eco_year, sumallPC1, sum_meanT, sum_GDD, GSLthreedayneg3C, sum_precip, moisturedeficit))
colnames(clim) <- c("Year", "PC1","Temp","GDD","GSL","Precip","MD")

# Add clim data to commtype summary dataframe 
anpp_clim <- inner_join(x = clim, y = anpp_type_year_commtype, by = c("Year" = "year"))

# Add clim data to all plot dataframe
sdlprod_clim <- left_join(x = clim, y = df3, by = c("Year" = "year"))

# Add all extended summer data to all plot dataframe
sdlprod_clim_tot <- left_join(x = extsum, y = df3, by = c("eco_year" = "year"))

# Add all extended summer data to summarized dataframe
sdlprod_clim_tot <- left_join(x = extsum, y = anpp_type_year_commtype, 
                              by = c("eco_year" = "year"))

# All Subsets Regression for WM
sdlprod_clim_tot_WM <- subset(sdlprod_clim_tot, veg_class == "WM")
library(bestglm)
env <- select(sdlprod_clim_tot_WM, c(sumallPC1, sumallPC2, sum_meanT, sum_GDD, fivedayrunning5C, fivedayrunning12C,GSLthreedayneg3C,sum_precip,sum_moisturedeficit,sum_PET,iceoff_GL4))
M <- cor(env)
corrplot(M, method = "number", type = "lower")
y <- sdlprod_clim_tot_WM$NPP
Xy <- as.data.frame(cbind(env,y))
bestmodel <- bestglm(Xy, IC = "AIC", RequireFullEnumerationQ = TRUE, TopModels=10)
bestmodel
bestmodel$BestModels
nullmodel <- lm(NPP ~ 1, data = sdlprod_clim_tot_WM)
bestmodel <- lm(NPP ~ GSLthreedayneg3C + sum_moisturedeficit, 
                data = sdlprod_clim_tot_WM)
summary(bestmodel)
plot(bestmodel)
AIC(nullmodel)
AIC(bestmodel)
AIC(bestmodel) - AIC(nullmodel)
anova(nullmodel, bestmodel)
# Partial R2
library(rsq)
rsq.partial(bestmodel, adj = TRUE)
rsq.partial(bestmodel, adj = FALSE)

# Do regressions on the big dataframe
# PC1
dm1 <- lm(NPP ~ PC1, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ PC1, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ PC1, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # Sig
sb1 <- lm(NPP ~ PC1, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ PC1, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Sig

# Temp
dm1 <- lm(NPP ~ Temp, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ Temp, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # Sig
sb1 <- lm(NPP ~ Temp, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Sig

# GDD
dm1 <- lm(NPP ~ GDD, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ GDD, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ GDD, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # Sig
sb1 <- lm(NPP ~ GDD, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ GDD, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Sig

# GSL
dm1 <- lm(NPP ~ GSL, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ GSL, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ GSL, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # Sig
sb1 <- lm(NPP ~ GSL, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ GSL, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Sig

# Precip
dm1 <- lm(NPP ~ Precip, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # R2 = 0.08, p < 0.001
mm1 <- lm(NPP ~ Precip, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # Sig
wm1 <- lm(NPP ~ Precip, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Precip, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # Sig
ff1 <- lm(NPP ~ Precip, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Marginal

# MD
dm1 <- lm(NPP ~ MD, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # R2 = 0.01, p = 0.04
mm1 <- lm(NPP ~ MD, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ MD, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ MD, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # Sig
ff1 <- lm(NPP ~ MD, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # NS

# Interactions
dm1 <- lm(NPP ~ Temp*GSL, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # NS
mm1 <- lm(NPP ~ Temp*GSL, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp*GSL, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp*GSL, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*GSL, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # NS

dm1 <- lm(NPP ~ Temp*Precip, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # Temp, Precip and Interaction Significant
mm1 <- lm(NPP ~ Temp*Precip, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp*Precip, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp*Precip, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*Precip, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # NS

dm1 <- lm(NPP ~ Temp*MD, data = subset(anpp_clim, veg_class == "DM"))
summary(dm1) # Significant interaction
mm1 <- lm(NPP ~ Temp*MD, data = subset(anpp_clim, veg_class == "MM"))
summary(mm1) # NS
wm1 <- lm(NPP ~ Temp*MD, data = subset(anpp_clim, veg_class == "WM"))
summary(wm1) # NS
sb1 <- lm(NPP ~ Temp*MD, data = subset(anpp_clim, veg_class == "SB"))
summary(sb1) # NS
ff1 <- lm(NPP ~ Temp*MD, data = subset(anpp_clim, veg_class == "FF"))
summary(ff1) # Temp
