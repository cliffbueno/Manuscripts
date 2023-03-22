# Plot MAG abundance across the 24 samples
library(plyr)
library(tidyverse)
library(reshape2)
library(car)
library(PMCMR)
meta <- read.csv("~/Desktop/Wetlands/GenomeSize.csv")
d <- read.delim("~/Desktop/Wetlands/data/Genomic_bins/bin_read_count.txt") %>%
  column_to_rownames(var = "I.OTU_ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample") %>%
  separate(Sample, into = c("junk", "ID"), sep = "MetaG") %>%
  select(-junk) %>%
  separate(ID, into = c("ID", "junk"), sep = "_MG") %>%
  select(-junk) %>%
  mutate(ID = str_replace(ID, "_R", "R")) %>%
  filter(ID %in% meta$ID) %>%
  select(ID, Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21) %>%
  left_join(., meta, by = "ID") %>%
  mutate(CPM = (Salt_Pond_MetaG_R1_A_D2_MG_DASTool_bins_concoct_out.21 * 1000000)/AssembledSize) %>%
  mutate(Site = gsub("R2A ", "R2A", Site)) %>%
  mutate(Site = factor(Site, levels = c("R1", "R2", "SF2", "R2A")))

# MAG Abundance by Site
stat_lab <- data.frame(x = levels(d$Site),
                       y = c(135, 45, 10, 10),
                       label = c("a", "ab", "bc", "c"))
png("~/Desktop/Methanolobus/MAG_Site.png", width = 5, height = 4, unit = "in", res = 300)
ggplot(d, aes(Site, CPM, colour = Restoration)) +
  geom_boxplot() +
  geom_jitter(size = 3, width = 0.2) +
  geom_text(data = stat_lab, aes(x = x, y = y, label = label),
            inherit.aes = F, show.legend = F, colour = "black", size = 5) +
  labs(x = "Site",
       y = "Counts per million",
       colour = "Wetland status") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(colour = "black", size = 0.25),
        axis.title = element_text(size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"))
dev.off()

# MAG abundance by methane and salinity
d_long <- melt(d, id.vars = c("Restoration", "CPM"), measure.vars = c("Methane", "Salinity"))
facet_names <- c("Methane" = "Methane flux (µmol/m^2/day)",
                 "Salinity" = "Salinity (ppt Cl)")
ggplot(d_long, aes(value, CPM, colour = Restoration)) +
  geom_point(size = 3) +
  geom_smooth(aes(value, CPM), se = F, inherit.aes = F, size = 0.25) +
  labs(x = NULL,
       y = "Counts per million assembled reads",
       colour = "Wetland status") +
  scale_colour_viridis_d() +
  facet_wrap(~ variable, ncol = 2, scales = "free_x", labeller = as_labeller(facet_names)) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(colour = "black", size = 0.25),
        axis.title = element_text(size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"))

# Methane and Salinity
ggplot(d, aes(Salinity, Methane, colour = Restoration)) +
  geom_point(size = 3) +
  geom_smooth(aes(Salinity, Methane), method = "lm", formula = y ~ x + I(x^2), 
              inherit.aes = F, size = 0.25, se = F) +
  labs(y = expression(paste("Methane flux (µmol/", m^{2}, "/day)")),
       x = "Salinity (Cl ppt)",
       colour = "Wetland status") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(colour = "black", size = 0.25),
        axis.title = element_text(size = 12, color = "black"), 
        axis.text = element_text(size = 10, color = "black"),
        axis.ticks = element_line(color = "black"))

# MAG abundance, Site - Methane - Salinity
d_long2 <- melt(d, id.vars = c("Restoration", "CPM"), measure.vars = c("Site","Methane", "Salinity")) %>%
  mutate(Site = value) %>%
  mutate(value = gsub("R1", "1", value)) %>%
  mutate(value = gsub("R2", "2", value)) %>%
  mutate(value = gsub("SF2", "3", value)) %>%
  mutate(value = gsub("2A", "4", value)) %>%
  mutate(value = as.numeric(value))
d_long2$variable <- factor(d_long2$variable,
                           levels = c("Site", "Methane", "Salinity"),
                           labels = c("(a)~Site",
                                      "(b)~Methane~flux~(µmol/m^{2}/day)",
                                      "(c)~Salinity~(ppt~chloride)"))
site_lab <- data.frame(x = c(1,2,3,4,0,500,1000,1500,50,100,150,200),
                       y = c(-13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13,-13),
                       variable = c("(a)~Site", 
                                    "(a)~Site", 
                                    "(a)~Site", 
                                    "(a)~Site",
                                    "(b)~Methane~flux~(µmol/m^{2}/day)", 
                                    "(b)~Methane~flux~(µmol/m^{2}/day)", 
                                    "(b)~Methane~flux~(µmol/m^{2}/day)", 
                                    "(b)~Methane~flux~(µmol/m^{2}/day)",
                                    "(c)~Salinity~(ppt~chloride)", 
                                    "(c)~Salinity~(ppt~chloride)", 
                                    "(c)~Salinity~(ppt~chloride)", 
                                    "(c)~Salinity~(ppt~chloride)"),
                       label = c("R1", "R2", "SF2", "R2A",
                                 "0", "500", "1000", "1500",
                                 "50", "100", "150", "200"))
stat_lab <- data.frame(x = c(1,2,3,4),
                       y = c(130, 130, 130, 130),
                       variable = c("(a)~Site", "(a)~Site", "(a)~Site", "(a)~Site"),
                       label = c("a", "ab", "bc", "c"))
site_lab$variable <- as.factor(site_lab$variable)
pdf("~/Desktop/Methanolobus/Manuscript/MAG_Abundance.pdf", width = 8, height = 3)
ggplot(d_long2, aes(value, CPM, colour = Restoration)) +
  geom_point(data = subset(d_long2, variable != "(a)~Site"),
             size = 3, alpha = 0.75) +
  geom_jitter(data = subset(d_long2, variable == "(a)~Site"),
              size = 3, alpha = 0.75, width = 0.2) +
  geom_smooth(data = subset(d_long2, variable != "(a)~Site"),
              aes(value, CPM), se = F, inherit.aes = F, size = 0.25) +
  geom_label(data = site_lab, aes(x = x, y = y, label = label),
            inherit.aes = F, show.legend = F, 
            colour = "black", fill = "white", size = 3.5, label.size = NA) +
  geom_text(data = stat_lab, aes(x = x, y = y, label = label),
            inherit.aes = F, show.legend = F, colour = "black", size = 3.5) +
  labs(x = NULL,
       y = "Counts per million",
       colour = "Wetland status") +
  scale_colour_viridis_d() +
  facet_wrap(~ variable, ncol = 3, scales = "free_x", labeller = label_parsed) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(colour = "black", size = 0.5),
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(size = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, colour = "black", hjust = 0),
        plot.margin = margin(0.1, 0.25, 1, 0.25, "lines")) +
  coord_cartesian(clip = "off",
                  ylim = c(0, 130))
dev.off()

# Stats
d$Site <- as.factor(d$Site)
leveneTest(CPM ~ Site, data = d)
shapiro.test(d$CPM)
m <- aov(CPM ~ Site, data = d)
shapiro.test(m$residuals)
summary(m)
TukeyHSD(m)

kruskal.test(CPM ~ Site, data = d)
posthoc.kruskal.nemenyi.test(CPM ~ Site, data = d)
