# ANI Analysis
# Make barplots and clustering diagrams for FastANI Kbase results

# Setup
library(plyr)
library(tidyverse)
library(readxl)
library(reshape2)
setwd("~/Desktop/Marivita/")
name_key <- read_xlsx("Marivita_FastANI.xlsx", sheet = 2)
d <- read_xlsx("Marivita_FastANI.xlsx", sheet = 1) %>%
  mutate_all(funs(str_replace(., name_key$ID[1], name_key$Name[1]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[2], name_key$Name[2]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[3], name_key$Name[3]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[4], name_key$Name[4]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[5], name_key$Name[5]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[6], name_key$Name[6]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[7], name_key$Name[7]))) %>%
  mutate_all(funs(str_replace(., "_assembly", "")))
d_wide <- dcast(d, QUERY ~ REFERENCE, value.var = "ANI_ESTIMATE") %>%
  column_to_rownames(var = "QUERY") %>%
  mutate_if(is.character, as.numeric) %>%
  as.matrix()
d_wide[upper.tri(d_wide)] <- NA
dist <- as.dist(d_wide)

# Cluster by ANI
x <- as.dendrogram(hclust(dist, method = "ward.D2"))
par(oma = c(1,1,1,1),
    mar = c(2,2,2,10))
plot(x, horiz = TRUE)

# Plot pairwise ANI with MAGs
MAGs_ANI <- d %>%
  filter(QUERY == "M. sp. SBSPR1" | QUERY == "M. sp. SBSPR2") %>%
  mutate(REFERENCE = as.factor(REFERENCE),
         ANI_ESTIMATE = as.numeric(ANI_ESTIMATE))

pdf("~/Desktop/Marivita/Manuscript/FigureS2.pdf", width = 7, height = 5)
ggplot(MAGs_ANI, aes(reorder(REFERENCE, ANI_ESTIMATE, median), ANI_ESTIMATE)) +
  geom_bar(stat = "identity") +
  labs(x = "Reference genome",
       y = "Average nucleotide identity (%)") +
  coord_flip(ylim = c(79, 85)) +
  facet_wrap(~ QUERY, ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black", family = "sans"), 
        axis.text.x = element_text(size = 10, color = "black", family = "sans"),
        axis.text.y = element_text(size = 10, color = "black", family = "sans", face = "italic"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face = "italic", size = 10))
dev.off()
