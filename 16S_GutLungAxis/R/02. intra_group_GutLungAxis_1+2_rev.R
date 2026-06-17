rm(list = ls())

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week1")
library("ggplot2")
library("ggpubr")
library(rstatix)
library(dplyr)

uw_raw <- read.csv("uw_dist.tsv", header = TRUE, sep = "\t")
w_raw  <- read.csv("w_dist.tsv",  header = TRUE, sep = "\t")

#####################################################################################
# (핵심) intra distance만: Group1 == Group2
uw <- uw_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Control","VNAM")) %>%
  transmute(Group = factor(Group1, levels = c("VNAM","Control")),
            Distance = Distance)

w <- w_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Control","VNAM")) %>%
  transmute(Group = factor(Group1, levels = c("VNAM","Control")),
            Distance = Distance)

#####################################################################################
## Unweighted UniFrac (INTRA)

p_df <- uw %>%
  wilcox_test(Distance ~ Group) %>%
  add_significance("p") %>%
  mutate(
    y.position = max(uw$Distance, na.rm = TRUE) * 1.05,
    g1   = as.numeric(factor(group1, levels = levels(uw$Group))),
    g2   = as.numeric(factor(group2, levels = levels(uw$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf("%s vs. %s : p = %.3g",
                      p_df$group1[1], p_df$group2[1], p_df$p[1])

intra_uw <- ggplot(uw, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position = "none",
        plot.title = element_text(size=24, hjust=0),
        axis.text.x = element_text(size=21),
        axis.text.y = element_text(size=24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c('royalblue','firebrick')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Unweighted UniFrac] Intra-group Distances") +
  stat_pvalue_manual(p_df,
                     label = NULL,
                     xmin  = "group1",
                     xmax  = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 0.3) +
  geom_text(data = p_df,
            inherit.aes = FALSE,
            aes(x = xpos * 1.05,
                y = y.position * 1.02,
                label = p.signif),
            angle = 270,
            vjust = 0.5,
            hjust = 0,
            size  = 9) +
  coord_flip(clip = "off") +
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1.1, vjust = 2.5, size = 6)

print(intra_uw)
ggsave("intra-group_rev/intra_uw.tiff", plot = intra_uw, width = 9, height = 7.5, dpi = 500)

#####################################################################################
## Weighted UniFrac (INTRA)

p_df <- w %>%
  wilcox_test(Distance ~ Group) %>%
  add_significance("p") %>%
  mutate(
    y.position = max(w$Distance, na.rm = TRUE) * 1.05,
    g1   = as.numeric(factor(group1, levels = levels(w$Group))),
    g2   = as.numeric(factor(group2, levels = levels(w$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf("%s vs. %s : p = %.3g",
                      p_df$group1[1], p_df$group2[1], p_df$p[1])

intra_w <- ggplot(w, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position = "none",
        plot.title = element_text(size=24, hjust=0),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c('royalblue','firebrick')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Weighted UniFrac] Intra-group Distances") +
  stat_pvalue_manual(p_df,
                     label = NULL,
                     xmin  = "group1",
                     xmax  = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 0.5) +
  geom_text(data = p_df,
            inherit.aes = FALSE,
            aes(x = xpos * 1.05,
                y = y.position * 1.02,
                label = p.signif),
            angle = 270,
            vjust = 0.5,
            hjust = 0,
            size  = 9) +
  coord_flip(clip = "off") +
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1.1, vjust = 2.5, size = 6)

print(intra_w)
ggsave("intra-group_rev/intra_w.tiff", plot = intra_w, width = 9, height = 7.5, dpi = 500)

#####################################################################################
#####################################################################################
rm(list = ls())

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week2")
library("ggplot2")
library("ggpubr")
library(rstatix)
library(dplyr)

bray_raw <- read.csv("bray_dist.tsv", header = TRUE, sep = "\t")
jaccard_raw  <- read.csv("jaccard_dist.tsv",  header = TRUE, sep = "\t")

#####################################################################################
# (핵심) intra distance만: Group1 == Group2
bray <- bray_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Control","VNAM")) %>%
  transmute(Group = factor(Group1, levels = c("VNAM","Control")),
            Distance = Distance)

jaccard <- jaccard_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Control","VNAM")) %>%
  transmute(Group = factor(Group1, levels = c("VNAM","Control")),
            Distance = Distance)

#####################################################################################
## Bray-Curtis (INTRA)

p_df <- bray %>%
  wilcox_test(Distance ~ Group) %>%
  add_significance("p") %>%
  mutate(
    y.position = max(bray$Distance, na.rm = TRUE) * 1.05,
    g1   = as.numeric(factor(group1, levels = levels(bray$Group))),
    g2   = as.numeric(factor(group2, levels = levels(bray$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf("%s vs. %s : p = %.3g",
                      p_df$group1[1], p_df$group2[1], p_df$p[1])

intra_bray <- ggplot(bray, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position = "none",
        plot.title = element_text(size=24, hjust=0),
        axis.text.x = element_text(size=21),
        axis.text.y = element_text(size=24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c('royalblue','firebrick')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Bray-Curtis Dissimilarity] Intra-group Distances") +
  stat_pvalue_manual(p_df,
                     label = NULL,
                     xmin  = "group1",
                     xmax  = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 0.3) +
  geom_text(data = p_df,
            inherit.aes = FALSE,
            aes(x = xpos * 1.05,
                y = y.position * 1.02,
                label = p.signif),
            angle = 270,
            vjust = 0.5,
            hjust = 0,
            size  = 9) +
  coord_flip(clip = "off") +
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1.1, vjust = 2.5, size = 6)

print(intra_bray)
ggsave("intra-group_rev/intra_bray.tiff", plot = intra_bray, width = 9, height = 7.5, dpi = 500)

#####################################################################################
## Jaccard (INTRA)

p_df <- jaccard %>%
  wilcox_test(Distance ~ Group) %>%
  add_significance("p") %>%
  mutate(
    y.position = max(jaccard$Distance, na.rm = TRUE) * 1.05,
    g1   = as.numeric(factor(group1, levels = levels(jaccard$Group))),
    g2   = as.numeric(factor(group2, levels = levels(jaccard$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf("%s vs. %s : p = %.3g",
                      p_df$group1[1], p_df$group2[1], p_df$p[1])

intra_jaccard <- ggplot(jaccard, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom ='errorbar', width = 0.4) +
  geom_boxplot() +
  theme_test() +
  theme(legend.position = "none",
        plot.title = element_text(size=24, hjust=0),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = c('royalblue','firebrick')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Jaccard Dissimilarity] Intra-group Distances") +
  stat_pvalue_manual(p_df,
                     label = NULL,
                     xmin  = "group1",
                     xmax  = "group2",
                     y.position = "y.position",
                     tip.length = 0.01,
                     size = 0.5) +
  geom_text(data = p_df,
            inherit.aes = FALSE,
            aes(x = xpos * 1.05,
                y = y.position * 1.02,
                label = p.signif),
            angle = 270,
            vjust = 0.5,
            hjust = 0,
            size  = 9) +
  coord_flip(clip = "off") +
  annotate("text",
           x = Inf, y = Inf,
           label = label_text,
           hjust = 1.1, vjust = 2.5, size = 6)

print(intra_jaccard)
ggsave("intra-group_rev/intra_jaccard.tiff", plot = intra_jaccard, width = 9, height = 7.5, dpi = 500)

#####################################################################################

