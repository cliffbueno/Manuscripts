# Calculate proteome median isoelectric point

library(readxl)
library(tidyverse)
setwd("~/Desktop/Methanolobus/")
se <- function(x) sd(x)/sqrt(length(x))

#### Isoelectric Point Calculator 2.0 ####
# Note only gets the first ~115 proteins (50k characters)!
# MAG: 5.55
d <- read_xlsx("IsoelectricPoint.xlsx")
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetBom: 5.56
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 2)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetPro: 5.63
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 3)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetPsyp: 5.86
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 4)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetPsyt: 5.92
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 5)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetSY01: 5.77 (note used average pI)
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 6)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('average pI', Col1)) %>%
  filter(., !grepl('includes all scales', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "Junk2", "Junk3", "pI"), sep = " ")
d$pI[69] <- 10.28
d$pI[81] <- 11.28
d <- d %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetTin: 5.38
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 7)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetVul: 5.36
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 8)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetVul1: 5.70
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 9)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

# MetZin: 5.37
d <- read_xlsx("IsoelectricPoint.xlsx", sheet = 10)
names(d)[1] <- "Col1"
d <- d %>%
  filter(., grepl('Isoelectric point of the protein is:', Col1)) %>%
  select(1) %>%
  separate(., Col1, into = c("Junk", "pI"), sep = ": ") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI)
hist(d$pI)

#### Galaxy ####
# Calculates for entire proteome! Use this.
setwd("~/Desktop/Methanolobus/IEP/")

# MAG
d <- read.delim("MAG.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
se(d$pI)
median(d$pI) # 5.00
hist(d$pI)

# MetBom
d <- read.delim("MetBom.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.09
hist(d$pI)

# MetPro
d <- read.delim("MetPro.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 4.92
hist(d$pI)

# MetPsyp
d <- read.delim("MetPsyp.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.82
hist(d$pI)

# MetPsyt
d <- read.delim("MetPsyt.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.49
hist(d$pI)

# MetSY1
d <- read.delim("MetSY1.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.40
hist(d$pI)

# MetTin
d <- read.delim("MetTin.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.02
hist(d$pI)

# MetVul
d <- read.delim("MetVul.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.02
hist(d$pI)

# MetVul1
d <- read.delim("MetVul1.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 5.12
hist(d$pI)

# MetZin
d <- read.delim("MetZin.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 4.97
hist(d$pI)

# MetvHol
d <- read.delim("MetvHol.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI))
mean(d$pI)
median(d$pI) # 4.97
hist(d$pI)

#### Graph ####
d1 <- read.delim("MAG.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ms. sp. SBSPR1A")
d2 <- read.delim("MetBom.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. bombayensis")
d3 <- read.delim("MetPro.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. profundi")
d4 <- read.delim("MetPsyp.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. psychrophilus")
d5 <- read.delim("MetPsyt.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. psychrotolerans")
d6 <- read.delim("MetSY1.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. sp SY-01")
d7 <- read.delim("MetTin.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. tindarius")
d8 <- read.delim("MetVul.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. vulcani B1d")
d9 <- read.delim("MetVul1.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. vulcani PL 12/M")
d10 <- read.delim("MetZin.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Ml. zinderi")
d11 <- read.delim("MetvHol.txt", header = F) %>%
  filter(., grepl('Isoelectric Point =', V1)) %>%
  separate(., V1, into = c("Junk", "pI"), sep = "=") %>%
  mutate("pI" = as.numeric(pI)) %>%
  mutate("Species" = "Mmv. hollandica")

d_long <- rbind(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11)

pdf("~/Desktop/Methanolobus/Manuscript/FigS8.pdf", width = 6, height = 5)
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
