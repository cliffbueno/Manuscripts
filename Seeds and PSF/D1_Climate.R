# Analyze Recent D1 Precip and Temperature to report in paper methods

# Clean up enviro, read in needed libraries
rm(list=ls())
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/cliffbueno/Functions/master/Summary.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
library(readxl)
library(tidyverse)
library(vegan)
library(naniar)
library(FSA)
options(stringsAsFactors = F)
theme_set(theme_bw())
na_vals <- c(" ", "", NA, NaN, "NA", "NaN", ".")

# -- FUNCTIONS -----
# set up functions to read in tabular datasets from EDI dynamically
# SCE + CTW code to determine most recent version of package ID and read in current dataset on EDI
# function to determine current version of data package on EDI
getCurrentVersion<-function(edi_id){
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

# D1 Precipitation and Temperature Data

# Infilled precipitation data for D1 chart recorder, 1952 - ongoing, daily
# Use 1952 - 2020
D1_precip <- getTabular(186) %>%
  separate(., date, sep = "-", into = c("year","month","day")) %>%
  filter(year > 1998) %>%
  filter(year < 2019) %>%
  group_by(year) %>%
  summarise(ann_precip = sum(precip)) %>%
  ungroup()
mean(D1_precip$ann_precip)

# Air temperature data for D1 chart recorder, 1952 - ongoing
# Use 1952 - 2020
D1_temp <- getTabular(412)

# Infilled air temperature data for D1 chart recorder, 1952 - 2018, daily
D1_temp <- getTabular(187) %>%
  separate(., date, sep = "-", into = c("year","month","day")) %>%
  filter(year > 1998) %>%
  group_by(year) %>%
  summarise(ann_meantemp = mean(mean_temp)) %>%
  ungroup()
mean(D1_temp$ann_meantemp)
