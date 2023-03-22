# ANI Analysis
# Make barplots and clustering diagrams for FastANI Kbase results

# Setup
library(plyr)
library(tidyverse)
library(readxl)
library(reshape2)
library(gplots)
library(scales)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(vegan)
library(RColorBrewer)
setwd("~/Desktop/Methanolobus/")
name_key <- read_xlsx("FastANI.xlsx", sheet = 4)
d <- read_xlsx("FastANI.xlsx", sheet = 3) %>%
  mutate_all(funs(str_replace(., name_key$ID[1], name_key$Name[1]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[2], name_key$Name[2]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[3], name_key$Name[3]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[4], name_key$Name[4]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[5], name_key$Name[5]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[6], name_key$Name[6]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[7], name_key$Name[7]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[8], name_key$Name[8]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[9], name_key$Name[9]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[10], name_key$Name[10]))) %>%
  mutate_all(funs(str_replace(., name_key$ID[11], name_key$Name[11]))) %>%
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

# Plot pairwise ANI with MAG
MAG_ANI <- d %>%
  filter(QUERY == "MAG") %>%
  mutate(REFERENCE = as.factor(REFERENCE),
         ANI_ESTIMATE = as.numeric(ANI_ESTIMATE))

pdf("~/Desktop/Methanolobus/Manuscript/FigS5.pdf", width = 7, height = 5)
ggplot(MAG_ANI, aes(reorder(REFERENCE, ANI_ESTIMATE, median), ANI_ESTIMATE)) +
  geom_bar(stat = "identity") +
  labs(x = "Reference genome",
       y = "Average nucleotide identity (%)") +
  coord_flip(ylim = c(76.5, 78)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12, color = "black", family = "sans"), 
        axis.text.x = element_text(size = 10, color = "black", family = "sans"),
        axis.text.y = element_text(size = 10, color = "black", family = "sans", face = "italic"),
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_blank())
dev.off()

# Extract pairwise distances to make heatmap like 16S
d$ANI <- round(as.numeric(d$ANI_ESTIMATE), digits = 2)
d1 <- select(d, QUERY, REFERENCE, ANI)
MAG <- subset(d1, QUERY == "MAG")
MetZin <- subset(d1, QUERY == "Ml. zinderi")
MetPsyt <- subset(d1, QUERY == "Ml. psychrotolerans")
MetVulPL <- subset(d1, QUERY == "Ml. vulcani PL 12/M")
MetVulB1d <- subset(d1, QUERY == "Ml. vulcani B1d")
MetPro <- subset(d1, QUERY == "Ml. profundi")
MetSY <- subset(d1, QUERY == "Ml. sp. SY-01")
MetPsyp <- subset(d1, QUERY == "Ml. psychrophilus")
MetTin <- subset(d1, QUERY == "Ml. tindarius")
MetBom <- subset(d1, QUERY == "Ml. bombayensis")
MetHol <- subset(d1, QUERY == "Mmv. hollandica")

a <- read_xlsx("GenomeInfo.xlsx", sheet = 6) %>%
  column_to_rownames(var = "Sequence") %>%
  as.matrix()
f <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
mat <- f(a)

# Pretty heatmap
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

phm1 <- pheatmap(mat,
                 legend = F,
                 scale = "none",
                 angle_col = 315,
                 fontsize = 6,
                 fontsize_row = 8,
                 fontsize_col = 8,
                 cluster_rows = T,
                 cluster_cols = T,
                 na_col = "white",
                 display_numbers = T)
save_pheatmap_pdf(phm1, "ANI_Heatmap.pdf")


#### 16S #####
# Update - also plot heatmap of 16S % identities
s <- read_xlsx("GenomeInfo.xlsx", sheet = 5) %>%
  column_to_rownames(var = "Sequence") %>%
  as.matrix()
f <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
mat <- f(s)
mat[1,1] <- 100

# Pretty heatmap
save_pheatmap_pdf <- function(x, filename, width = 8, height = 8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

phm1 <- pheatmap(mat,
                 legend = F,
                 scale = "none",
                 angle_col = 315,
                 fontsize = 6,
                 fontsize_row = 8,
                 fontsize_col = 8,
                 cluster_rows = T,
                 cluster_cols = T,
                 na_col = "white",
                 display_numbers = T)
setwd("~/Desktop/Methanolobus/Manuscript/")
save_pheatmap_pdf(phm1, "FigS4.pdf")

