# Calculate proteome median isoelectric point

library(readxl)
library(tidyverse)
setwd("~/Desktop/Marivita/MarivitaIEP/")
se <- function(x) sd(x)/sqrt(length(x))

#### Galaxy ####
# Calculates for entire proteome! Use this.

# MarCry
d <- read.delim("MarCry.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
se(d$pI)
median(d$pI)
hist(d$pI)

# MarGeo
d <- read.delim("MarGeo.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MarHal
d <- read.delim("MarHal.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# Marivita 10_192
d <- read.delim("Marivita10_192.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# Marivita 28_82
d <- read.delim("Marivita28_82.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MarLac
d <- read.delim("MarLac.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MarLZ
d <- read.delim("MarLZ.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)


#### Graph ####
d1 <- read.delim("MarCry.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M.cryptomonadis MP20-4")
d2 <- read.delim("MarGeo.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. geojedonensis")
d3 <- read.delim("MarHal.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. hallyeonensis")
d4 <- read.delim("Marivita10_192.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. sp. SBSPR2")
d5 <- read.delim("Marivita28_82.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. sp. SBSPR1")
d6 <- read.delim("MarLac.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. lacus")
d7 <- read.delim("MarLZ.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "M. cryptomonadis LZ-15-2")

d_long <- rbind(d1, d2, d3, d4, d5, d6, d7)

pdf("~/Desktop/Marivita/Manuscript/FigureS7.pdf", width = 7, height = 4)
ggplot() +
  geom_histogram(data = d_long, aes(x = pI)) +
  labs(x = "Isoelectric point",
         y = "Number of proteins") +
  facet_wrap(~ Species, ncol = 4) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black"), 
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.text = element_text(size = 8, color = "black", face = "italic"))
dev.off()
