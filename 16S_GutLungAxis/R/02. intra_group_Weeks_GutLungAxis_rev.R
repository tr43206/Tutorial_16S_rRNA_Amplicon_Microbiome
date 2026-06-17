rm(list = ls())

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_VNAM")
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(mutoss)   # BKY
#install.packages("mutoss", repos = "https://cloud.r-project.org")

uw_raw <- read.csv("uw_dist.tsv", header=TRUE, sep="\t")
w_raw  <- read.csv("w_dist.tsv",  header=TRUE, sep="\t")

#####################################################################################
# (핵심) intra distance만: Group1 == Group2
uw <- uw_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Week0","Week1","Week2")) %>%
  transmute(Group = factor(Group1, levels = c("Week2","Week1","Week0")),
            Distance = Distance)

w <- w_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Week0","Week1","Week2")) %>%
  transmute(Group = factor(Group1, levels = c("Week2","Week1","Week0")),
            Distance = Distance)

# 비교쌍
comp <- list(c("Week0","Week1"), c("Week0","Week2"), c("Week1","Week2"))

# BKY 보정 함수
p_adjust_BKY <- function(p, max_iter = 60){
  p <- as.numeric(p)
  n <- length(p)
  out <- rep(NA_real_, n)
  
  ok <- which(!is.na(p))
  if(length(ok) == 0) return(out)
  
  # ok subset 정렬
  ord <- order(p[ok])
  p_sorted <- p[ok][ord]
  m <- length(p_sorted)
  
  k2_for_q <- function(q){
    if(q <= 0) return(0L)
    q1 <- q / (1 + q)                  # stage1 level
    thr1 <- (seq_len(m) * q1 / m)      # BH thresholds at q1
    k1 <- if(any(p_sorted <= thr1)) max(which(p_sorted <= thr1)) else 0L
    m0 <- max(1L, m - k1)              # estimated # true nulls
    thr2 <- (seq_len(m) * q1 / m0)     # stage2 thresholds
    k2 <- if(any(p_sorted <= thr2)) max(which(p_sorted <= thr2)) else 0L
    k2
  }
  
  k2_at_1 <- k2_for_q(1)
  adj_sorted <- rep(1, m)
  
  for(i in seq_len(m)){
    if(i > k2_at_1){ adj_sorted[i] <- 1; next }
    lo <- 0; hi <- 1
    for(it in seq_len(max_iter)){
      mid <- (lo + hi) / 2
      if(i <= k2_for_q(mid)) hi <- mid else lo <- mid
    }
    adj_sorted[i] <- hi
  }
  
  # rank별 monotone 보정
  adj_sorted <- cummax(adj_sorted)
  adj_ok <- numeric(m); adj_ok[ord] <- adj_sorted
  
  out[ok] <- pmin(adj_ok, 1)
  out
}


# p값 뽑기(라벨용) - BKY 보정 p 사용
get_p <- function(df, a, b, col = "p.adj"){
  out <- df[[col]][(df$group1==a & df$group2==b) | (df$group1==b & df$group2==a)]
  ifelse(length(out)==0, NA, out[1])
}

#####################################################################################
## Unweighted UniFrac (INTRA) : Kruskal-Wallis + BKY post-hoc

kw_uw   <- kruskal_test(uw, Distance ~ Group)
kw_p_uw <- kw_uw$p

y_max <- max(uw$Distance, na.rm=TRUE)

p_df <- uw %>%
  pairwise_wilcox_test(Distance ~ Group, comparisons = comp, p.adjust.method = "none") %>%
  mutate(p.adj = p_adjust_BKY(p)) %>%         # <- BKY 보정 p
  add_significance("p.adj") %>%             # p.adj.signif 생성
  mutate(
    p.signif   = p.adj.signif,              # <- 별표는 BKY 기준
    y.position = y_max * c(1.05, 1.15, 1.25),
    g1   = as.numeric(factor(group1, levels = levels(uw$Group))),
    g2   = as.numeric(factor(group2, levels = levels(uw$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf(
  "Week0 vs. Week1 : q = %.3g\nWeek0 vs. Week2 : q = %.3g\nWeek1 vs. Week2 : q = %.3g",
  get_p(p_df,"Week0","Week1","p.adj"),
  get_p(p_df,"Week0","Week2","p.adj"),
  get_p(p_df,"Week1","Week2","p.adj")
)

intra_uw <- ggplot(uw, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom='errorbar', width=0.4) +
  geom_boxplot() +
  theme_test() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=22, hjust=0),
    axis.text.x = element_text(size=21),
    axis.text.y = element_text(size=24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 35, 10, 10, "pt")
  ) +
  scale_fill_manual(values=c('#008080','#0D98BA','gray')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Unweighted UniFrac] Intra-group Distances") +
  stat_pvalue_manual(
    p_df, label=NULL,
    xmin="group1", xmax="group2",
    y.position="y.position",
    tip.length=0.01, size=0.3
  ) +
  geom_text(
    data=p_df, inherit.aes=FALSE,
    aes(x=xpos*1.025, y=y.position*1.03, label=p.adj.signif),
    angle=270, vjust=0.5, hjust=0, size=8
  ) +
  coord_flip(clip="off") +
  annotation_custom(
    grid::textGrob(
      label_text,
      x = unit(0.6, "npc"),  # 왼쪽/오른쪽 위치 (0~1)
      y = unit(0.98, "npc"),  # 위/아래 위치 (0~1)
      just = c("left", "top"),# <- 핵심: 왼쪽 정렬 + 위 기준
      gp = grid::gpar(fontsize = 14)
    )
  )

print(intra_uw)
ggsave("intra-group_rev/weeks_intra_uw.tiff", plot=intra_uw, width=9, height=7.5, dpi=500)

#####################################################################################
## Weighted UniFrac (INTRA) : Kruskal-Wallis + BKY post-hoc

kw_w   <- kruskal_test(w, Distance ~ Group)
kw_p_w <- kw_w$p

y_max <- max(w$Distance, na.rm=TRUE)

p_df <- w %>%
  pairwise_wilcox_test(Distance ~ Group, comparisons = comp, p.adjust.method = "none") %>%
  mutate(p.adj = p_adjust_BKY(p)) %>%         # <- BKY 보정 p
  add_significance("p.adj") %>%             # p.adj.signif 생성
  mutate(
    p.signif   = p.adj.signif,              # <- 별표는 BKY 기준
    y.position = y_max * c(1.05, 1.15, 1.25),
    g1   = as.numeric(factor(group1, levels = levels(w$Group))),
    g2   = as.numeric(factor(group2, levels = levels(w$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf(
  "Week0 vs. Week1 : q = %.3g\nWeek0 vs. Week2 : q = %.3g\nWeek1 vs. Week2 : q = %.3g",
  get_p(p_df,"Week0","Week1","p.adj"),
  get_p(p_df,"Week0","Week2","p.adj"),
  get_p(p_df,"Week1","Week2","p.adj")
)

intra_w <- ggplot(w, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom='errorbar', width=0.4) +
  geom_boxplot() +
  theme_test() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=22, hjust=0),
    axis.text.x = element_text(size=21),
    axis.text.y = element_text(size=24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 35, 10, 10, "pt")
  ) +
  scale_fill_manual(values=c('#008080','#0D98BA','gray')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Weighted UniFrac] Intra-group Distances") +
  stat_pvalue_manual(
    p_df, label=NULL,
    xmin="group1", xmax="group2",
    y.position="y.position",
    tip.length=0.01, size=0.3
  ) +
  geom_text(
    data=p_df, inherit.aes=FALSE,
    aes(x=xpos*1.025, y=y.position*1.03, label=p.adj.signif),
    angle=270, vjust=0.5, hjust=0, size=8
  ) +
  coord_flip(clip="off") +
  annotation_custom(
    grid::textGrob(
      label_text,
      x = unit(0.6, "npc"),  # 왼쪽/오른쪽 위치 (0~1)
      y = unit(0.98, "npc"),  # 위/아래 위치 (0~1)
      just = c("left", "top"),# <- 핵심: 왼쪽 정렬 + 위 기준
      gp = grid::gpar(fontsize = 14)
    )
  )

print(intra_w)
ggsave("intra-group_rev/weeks_intra_w.tiff", plot=intra_w, width=9, height=7.5, dpi=500)

##########################################################################################################################################################################
rm(list = ls())

bray_raw <- read.csv("bray_dist.tsv", header=TRUE, sep="\t")
jaccard_raw  <- read.csv("jaccard_dist.tsv",  header=TRUE, sep="\t")

#####################################################################################
# (핵심) intra distance만: Group1 == Group2
bray <- bray_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Week0","Week1","Week2")) %>%
  transmute(Group = factor(Group1, levels = c("Week2","Week1","Week0")),
            Distance = Distance)

jaccard <- jaccard_raw %>%
  filter(Group1 == Group2, Group1 %in% c("Week0","Week1","Week2")) %>%
  transmute(Group = factor(Group1, levels = c("Week2","Week1","Week0")),
            Distance = Distance)

# 비교쌍(원하는 순서)
comp <- list(c("Week0","Week1"), c("Week0","Week2"), c("Week1","Week2"))

# BKY 보정 함수
p_adjust_BKY <- function(p, max_iter = 60){
  p <- as.numeric(p)
  n <- length(p)
  out <- rep(NA_real_, n)
  
  ok <- which(!is.na(p))
  if(length(ok) == 0) return(out)
  
  # ok subset 정렬
  ord <- order(p[ok])
  p_sorted <- p[ok][ord]
  m <- length(p_sorted)
  
  k2_for_q <- function(q){
    if(q <= 0) return(0L)
    q1 <- q / (1 + q)                  # stage1 level
    thr1 <- (seq_len(m) * q1 / m)      # BH thresholds at q1
    k1 <- if(any(p_sorted <= thr1)) max(which(p_sorted <= thr1)) else 0L
    m0 <- max(1L, m - k1)              # estimated # true nulls
    thr2 <- (seq_len(m) * q1 / m0)     # stage2 thresholds
    k2 <- if(any(p_sorted <= thr2)) max(which(p_sorted <= thr2)) else 0L
    k2
  }
  
  k2_at_1 <- k2_for_q(1)
  adj_sorted <- rep(1, m)
  
  for(i in seq_len(m)){
    if(i > k2_at_1){ adj_sorted[i] <- 1; next }
    lo <- 0; hi <- 1
    for(it in seq_len(max_iter)){
      mid <- (lo + hi) / 2
      if(i <= k2_for_q(mid)) hi <- mid else lo <- mid
    }
    adj_sorted[i] <- hi
  }
  
  # rank별 monotone 보정
  adj_sorted <- cummax(adj_sorted)
  adj_ok <- numeric(m); adj_ok[ord] <- adj_sorted
  
  out[ok] <- pmin(adj_ok, 1)
  out
}


# p값 뽑기(라벨용) - BKY 보정 p 사용
get_p <- function(df, a, b, col = "p.adj"){
  out <- df[[col]][(df$group1==a & df$group2==b) | (df$group1==b & df$group2==a)]
  ifelse(length(out)==0, NA, out[1])
}

#####################################################################################
## Bray-Curtis (INTRA) : Kruskal-Wallis + BKY post-hoc

kw_bray   <- kruskal_test(bray, Distance ~ Group)
kw_p_bray <- kw_bray$p

y_max <- max(bray$Distance, na.rm=TRUE)

p_df <- bray %>%
  pairwise_wilcox_test(Distance ~ Group, comparisons = comp, p.adjust.method = "none") %>%
  mutate(p.adj = p_adjust_BKY(p)) %>%         # <- BKY 보정 p
  add_significance("p.adj") %>%             # p.adj.signif 생성
  mutate(
    p.signif   = p.adj.signif,              # <- 별표는 BKY 기준
    y.position = y_max * c(1.05, 1.15, 1.25),
    g1   = as.numeric(factor(group1, levels = levels(bray$Group))),
    g2   = as.numeric(factor(group2, levels = levels(bray$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf(
  "Week0 vs. Week1 : q = %.3g\nWeek0 vs. Week2 : q = %.3g\nWeek1 vs. Week2 : q = %.3g",
  get_p(p_df,"Week0","Week1","p.adj"),
  get_p(p_df,"Week0","Week2","p.adj"),
  get_p(p_df,"Week1","Week2","p.adj")
)

intra_bray <- ggplot(bray, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom='errorbar', width=0.4) +
  geom_boxplot() +
  theme_test() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=22, hjust=0),
    axis.text.x = element_text(size=21),
    axis.text.y = element_text(size=24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 35, 10, 10, "pt")
  ) +
  scale_fill_manual(values=c('#008080','#0D98BA','gray')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Bray-Curtis Dissimilarity] Intra-group Distances") +
  stat_pvalue_manual(
    p_df, label=NULL,
    xmin="group1", xmax="group2",
    y.position="y.position",
    tip.length=0.01, size=0.3
  ) +
  geom_text(
    data=p_df, inherit.aes=FALSE,
    aes(x=xpos*1.025, y=y.position*1.03, label=p.adj.signif),
    angle=270, vjust=0.5, hjust=0, size=8
  ) +
  coord_flip(clip="off") +
  annotation_custom(
    grid::textGrob(
      label_text,
      x = unit(0.6, "npc"),  # 왼쪽/오른쪽 위치 (0~1)
      y = unit(0.98, "npc"),  # 위/아래 위치 (0~1)
      just = c("left", "top"),# <- 핵심: 왼쪽 정렬 + 위 기준
      gp = grid::gpar(fontsize = 14)
    )
  )

print(intra_bray)
ggsave("intra-group_rev/weeks_intra_bray.tiff", plot=intra_bray, width=9, height=7.5, dpi=500)

#####################################################################################
## Jaccard Dissimilarity (INTRA) : Kruskal-Wallis + BKY post-hoc

kw_jaccard   <- kruskal_test(jaccard, Distance ~ Group)
kw_p_jaccard <- kw_jaccard$p

y_max <- max(jaccard$Distance, na.rm=TRUE)

p_df <- jaccard %>%
  pairwise_wilcox_test(Distance ~ Group, comparisons = comp, p.adjust.method = "none") %>%
  mutate(p.adj = p_adjust_BKY(p)) %>%         # <- BKY 보정 p
  add_significance("p.adj") %>%             # p.adj.signif 생성
  mutate(
    p.signif   = p.adj.signif,              # <- 별표는 BKY 기준
    y.position = y_max * c(1.05, 1.15, 1.25),
    g1   = as.numeric(factor(group1, levels = levels(jaccard$Group))),
    g2   = as.numeric(factor(group2, levels = levels(jaccard$Group))),
    xpos = (g1 + g2) / 2
  )

label_text <- sprintf(
  "Week0 vs. Week1 : q = %.3g\nWeek0 vs. Week2 : q = %.3g\nWeek1 vs. Week2 : q = %.3g",
  get_p(p_df,"Week0","Week1","p.adj"),
  get_p(p_df,"Week0","Week2","p.adj"),
  get_p(p_df,"Week1","Week2","p.adj")
)

intra_jaccard <- ggplot(jaccard, aes(x = Group, y = Distance, fill = Group)) +
  stat_boxplot(geom='errorbar', width=0.4) +
  geom_boxplot() +
  theme_test() +
  theme(
    legend.position = "none",
    plot.title = element_text(size=22, hjust=0),
    axis.text.x = element_text(size=21),
    axis.text.y = element_text(size=24),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(10, 35, 10, 10, "pt")
  ) +
  scale_fill_manual(values=c('#008080','#0D98BA','gray')) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "[Jaccard Dissimilarity] Intra-group Distances") +
  stat_pvalue_manual(
    p_df, label=NULL,
    xmin="group1", xmax="group2",
    y.position="y.position",
    tip.length=0.01, size=0.3
  ) +
  geom_text(
    data=p_df, inherit.aes=FALSE,
    aes(x=xpos*1.025, y=y.position*1.03, label=p.adj.signif),
    angle=270, vjust=0.5, hjust=0, size=8
  ) +
  coord_flip(clip="off") +
  annotation_custom(
    grid::textGrob(
      label_text,
      x = unit(0.6, "npc"),  # 왼쪽/오른쪽 위치 (0~1)
      y = unit(0.98, "npc"),  # 위/아래 위치 (0~1)
      just = c("left", "top"),# <- 핵심: 왼쪽 정렬 + 위 기준
      gp = grid::gpar(fontsize = 14)
    )
  )

print(intra_jaccard)
ggsave("intra-group_rev/weeks_intra_jaccard.tiff", plot=intra_jaccard, width=9, height=7.5, dpi=500)

#####################################################################################
