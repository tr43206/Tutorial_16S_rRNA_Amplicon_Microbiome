############################################################
# Genus-level microbiome ordination for your data
# Input:
#   1) rel_gg2_genus_table(1).tsv
#   2) maaslin3_RUN45_Metadata_GutLungAxis(7).txt
#
# Analyses:
#   PCA
#   RDA
#   partial RDA
#   DCA
#   CCA
#   Bray-Curtis PCoA
#   PERMANOVA
#   dbRDA / capscale
############################################################

rm(list = ls())

############################################################
# 0. Packages
############################################################

pkg <- c("tidyverse", "vegan")

for (p in pkg) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(tidyverse)
library(vegan)

setwd('C:/Users/tr432/Downloads/filter_GutLungAxis_Week2')

############################################################
# 1. Input files
############################################################

taxa_file <- 'rel_gg2_genus_table.tsv'
meta_file <- 'C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt'

outdir <- 'ordination_results'
dir.create(outdir, showWarnings = FALSE)

############################################################
# 2. Read genus table
############################################################

# rel_gg2_genus_table은 첫 줄이 "# Constructed from biom file"이므로 skip = 1 필요
taxa_raw <- read.delim(
  taxa_file,
  sep = "\t",
  header = TRUE,
  skip = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

# 첫 번째 열 이름 정리
colnames(taxa_raw)[1] <- "Taxon"

# taxonomy label 정리 함수
clean_taxon_name <- function(x) {
  sapply(strsplit(x, ";"), function(v) {
    v <- gsub("^[a-z]__", "", v)
    v <- gsub("^__", "", v)
    v <- v[!is.na(v)]
    v <- v[v != "" & v != "__"]
    
    if (length(v) == 0) {
      return("Unclassified")
    } else {
      return(tail(v, 1))
    }
  }, USE.NAMES = FALSE)
}

taxon_full <- taxa_raw$Taxon
taxon_label <- make.unique(clean_taxon_name(taxon_full))

taxa_df <- taxa_raw[, -1, drop = FALSE]
rownames(taxa_df) <- taxon_label

# 숫자형 변환
taxa_df <- as.data.frame(
  lapply(taxa_df, function(x) as.numeric(as.character(x))),
  check.names = FALSE
)

rownames(taxa_df) <- taxon_label

# taxa x sample matrix
taxa_taxon_sample <- as.matrix(taxa_df)

# sample x taxa matrix
taxa_mat <- t(taxa_taxon_sample)

############################################################
# 3. Read metadata
############################################################

meta <- read.delim(
  meta_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "N/A"),
  comment.char = ""
)

# #SampleID 이름 정리
colnames(meta)[colnames(meta) == "#SampleID"] <- "SampleID"

############################################################
# 4. Match samples between genus table and metadata
############################################################

common_samples <- intersect(rownames(taxa_mat), meta$SampleID)

taxa_mat <- taxa_mat[common_samples, , drop = FALSE]

meta_sub <- meta %>%
  filter(SampleID %in% common_samples) %>%
  arrange(match(SampleID, common_samples))

stopifnot(all(meta_sub$SampleID == rownames(taxa_mat)))

############################################################
# 5. Metadata formatting
############################################################

meta_sub$Group <- factor(meta_sub$Group, levels = c("Control", "VNAM"))
meta_sub$SamplingWeek <- factor(meta_sub$SamplingWeek)

# R formula에서 쓰기 편하게 변수명 변경
meta_sub$InitialBW <- as.numeric(meta_sub$`Initial_b.w`)
meta_sub$SeqDepth <- as.numeric(meta_sub$SequencingDepth)
meta_sub$MGX_Reads_num <- as.numeric(meta_sub$MGX_Reads)

# 현재 genus table이 어떤 week인지 확인
cat("\n===== SamplingWeek in matched data =====\n")
print(table(meta_sub$SamplingWeek, useNA = "ifany"))

cat("\n===== Group in matched data =====\n")
print(table(meta_sub$Group, useNA = "ifany"))

############################################################
# 6. Data check
############################################################

cat("\n===== Raw matrix check =====\n")
cat("taxa_mat dimension: ", dim(taxa_mat), "\n")
cat("metadata dimension: ", dim(meta_sub), "\n")
cat("Sample matching: ", all(meta_sub$SampleID == rownames(taxa_mat)), "\n")

cat("\nSample-wise total abundance before normalization:\n")
print(round(rowSums(taxa_mat), 6))

cat("\nZero proportion before filtering:\n")
print(mean(taxa_mat == 0))

############################################################
# 7. Normalize relative abundance
############################################################

# 이미 relative abundance지만, sample sum을 1로 재보정
taxa_rel <- taxa_mat / rowSums(taxa_mat)
taxa_rel[is.na(taxa_rel)] <- 0

cat("\nSample-wise total abundance after normalization:\n")
print(round(rowSums(taxa_rel), 6))

############################################################
# 8. Taxa filtering
############################################################

# 전체 sample 중 20% 이상에서 검출된 genus만 사용
prev_cutoff <- 0.05

keep_taxa <- colMeans(taxa_rel > 0) >= prev_cutoff &
  colSums(taxa_rel) > 0

taxa_rel_filt <- taxa_rel[, keep_taxa, drop = FALSE]

cat("\n===== Filtering result =====\n")
cat("Before filtering: ", ncol(taxa_rel), " taxa\n")
cat("After filtering:  ", ncol(taxa_rel_filt), " taxa\n")
cat("Zero proportion after filtering: ", mean(taxa_rel_filt == 0), "\n")

# taxonomy map 저장
taxon_map <- data.frame(
  TaxonLabel = taxon_label,
  FullTaxonomy = taxon_full,
  stringsAsFactors = FALSE
)

write.csv(
  taxon_map,
  file = file.path(outdir, "taxon_label_map.csv"),
  row.names = FALSE
)

############################################################
# 9. Hellinger transformation
############################################################

# RDA/PCA용 변환
taxa_hell <- decostand(taxa_rel_filt, method = "hellinger")

############################################################
# 10. PCA
############################################################

pca_res <- prcomp(taxa_hell, center = TRUE, scale. = FALSE)

pca_var <- summary(pca_res)$importance[2, ] * 100

pca_df <- data.frame(
  SampleID = rownames(taxa_hell),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Group = meta_sub$Group,
  SamplingWeek = meta_sub$SamplingWeek,
  InitialBW = meta_sub$InitialBW
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Control" = '#B22222', "VNAM" = '#4169E1')) +
  labs(
    title = "PCA of genus profiles",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_classic(base_size = 14)

if (min(table(pca_df$Group)) >= 3) {
  p_pca <- p_pca + stat_ellipse(type = "t", linewidth = 0.7)
}

print(p_pca)

ggsave(
  filename = file.path(outdir, "PCA_genus_Hellinger.png"),
  plot = p_pca,
  width = 5,
  height = 4
)

write.csv(
  pca_df,
  file = file.path(outdir, "PCA_scores.csv"),
  row.names = FALSE
)

cat("\n===== PCA variance explained =====\n")
print(round(pca_var[1:10], 3))

############################################################
# 11. PCA loading / eigenvalue check
############################################################

E <- cov(taxa_hell)
evd <- eigen(E)

cat("\n===== PCA eigenvalues top 10 =====\n")
print(evd$values[1:min(10, length(evd$values))])

cat("\n===== Covariance matrix first 5 x 5 =====\n")
print(E[1:min(5, nrow(E)), 1:min(5, ncol(E))])

############################################################
# 12. RDA: Group + Initial body weight
############################################################

rda_res <- rda(
  taxa_hell ~ InitialBW + Group,
  data = meta_sub
)

cat("\n===== RDA result: taxa_hell ~ InitialBW + Group =====\n")
print(rda_res)

cat("\n===== RDA adjusted R2 =====\n")
print(RsquareAdj(rda_res))

cat("\n===== RDA ANOVA: full model =====\n")
print(anova(rda_res, permutations = 999))

cat("\n===== RDA ANOVA: by term, sequential =====\n")
print(anova(rda_res, by = "term", permutations = 999))

cat("\n===== RDA ANOVA: by margin =====\n")
print(anova(rda_res, by = "margin", permutations = 999))

cat("\n===== RDA ANOVA: by axis =====\n")
print(anova(rda_res, by = "axis", permutations = 999))


############################################################
# RDA plot with group separation + sample labels
# Output: PNG
############################################################

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

# 점수 추출
site_sc <- scores(rda_res, display = "sites", choices = 1:2, scaling = 2)
sp_sc   <- scores(rda_res, display = "species", choices = 1:2, scaling = 2)
bp_sc   <- scores(rda_res, display = "bp", choices = 1:2, scaling = 2)

# sample 점수
site_df <- data.frame(
  SampleID = rownames(site_sc),
  RDA1 = site_sc[, 1],
  RDA2 = site_sc[, 2],
  Group = meta_sub$Group
)

# species 점수
sp_df <- data.frame(
  Taxon = rownames(sp_sc),
  RDA1 = sp_sc[, 1],
  RDA2 = sp_sc[, 2]
)

# 설명변수 점수
bp_df <- data.frame(
  Variable = rownames(bp_sc),
  RDA1 = bp_sc[, 1],
  RDA2 = bp_sc[, 2]
)

# 변수명 정리
bp_df$Variable <- dplyr::recode(
  bp_df$Variable,
  "InitialBW" = "InitialBW",
  "GroupVNAM" = "Group"
)

# genus arrow 길이 기준 상위 일부만 표시
top_n_taxa <- 8

sp_df <- sp_df %>%
  mutate(arrow_len = sqrt(RDA1^2 + RDA2^2)) %>%
  arrange(desc(arrow_len)) %>%
  slice_head(n = top_n_taxa)

# 축 설명력
rda_var <- summary(rda_res)$cont$importance[2, 1:2] * 100

# 화살표 스케일 조정
site_range <- max(abs(c(site_df$RDA1, site_df$RDA2)), na.rm = TRUE)
sp_range   <- max(abs(c(sp_df$RDA1, sp_df$RDA2)), na.rm = TRUE)
bp_range   <- max(abs(c(bp_df$RDA1, bp_df$RDA2)), na.rm = TRUE)

sp_mul <- site_range / sp_range * 0.85
bp_mul <- site_range / bp_range * 0.75

sp_df <- sp_df %>%
  mutate(
    RDA1_plot = RDA1 * sp_mul,
    RDA2_plot = RDA2 * sp_mul
  )

bp_df <- bp_df %>%
  mutate(
    RDA1_plot = RDA1 * bp_mul,
    RDA2_plot = RDA2 * bp_mul
  )

# plot 범위
x_lim <- range(c(site_df$RDA1, sp_df$RDA1_plot, bp_df$RDA1_plot), na.rm = TRUE)
y_lim <- range(c(site_df$RDA2, sp_df$RDA2_plot, bp_df$RDA2_plot), na.rm = TRUE)

x_pad <- diff(x_lim) * 0.12
y_pad <- diff(y_lim) * 0.12

# 그룹 색/모양
group_cols <- c("Control" = "red", "VNAM" = "blue")
group_shapes <- c("Control" = 16, "VNAM" = 17)

# plot
p_rda_group_label <- ggplot() +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.5) +
  
  # explanatory variable arrows
  geom_segment(
    data = bp_df,
    aes(x = 0, y = 0, xend = RDA1_plot, yend = RDA2_plot),
    arrow = arrow(length = unit(0.22, "cm")),
    color = "red",
    linewidth = 0.8
  ) +
  geom_text_repel(
    data = bp_df,
    aes(x = RDA1_plot, y = RDA2_plot, label = Variable),
    color = "red",
    size = 4,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  # sample points
  geom_point(
    data = site_df,
    aes(x = RDA1, y = RDA2, color = Group, shape = Group),
    size = 2.8
  ) +
  
  # sample labels
  geom_text_repel(
    data = site_df,
    aes(x = RDA1, y = RDA2, label = SampleID, color = Group),
    size = 3.2,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = "gray70",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  # genus arrows
  geom_segment(
    data = sp_df,
    aes(x = 0, y = 0, xend = RDA1_plot, yend = RDA2_plot),
    arrow = arrow(length = unit(0.22, "cm")),
    color = "darkgreen",
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = sp_df,
    aes(x = RDA1_plot, y = RDA2_plot, label = Taxon),
    color = "darkgreen",
    size = 3.6,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  
  scale_color_manual(values = group_cols) +
  scale_shape_manual(values = group_shapes) +
  
  coord_cartesian(
    xlim = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
    ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)
  ) +
  labs(
    title = "RDA on Genus level",
    x = paste0("RDA1(", round(rda_var[1], 2), "%)"),
    y = paste0("RDA2(", round(rda_var[2], 2), "%)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.title = element_blank(),
    legend.position = c(0.90, 0.92),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

print(p_rda_group_label)

ggsave(
  filename = file.path(outdir, "RDA_on_Genus_level_grouped_labeled.png"),
  plot = p_rda_group_label,
  width = 7.5,
  height = 5.5,
  dpi = 500
)

############################################################
# 13. Partial RDA: Group effect after conditioning InitialBW
############################################################

rda_partial <- rda(
  taxa_hell ~ Group + Condition(InitialBW),
  data = meta_sub
)

cat("\n===== Partial RDA result: taxa_hell ~ Group + Condition(InitialBW) =====\n")
print(rda_partial)

cat("\n===== Partial RDA adjusted R2 =====\n")
print(RsquareAdj(rda_partial))

cat("\n===== Partial RDA ANOVA: full model =====\n")
print(anova(rda_partial, permutations = 999))

cat("\n===== Partial RDA ANOVA: by term =====\n")
print(anova(rda_partial, by = "term", permutations = 999))

cat("\n===== Partial RDA ANOVA: by axis =====\n")
print(anova(rda_partial, by = "axis", permutations = 999))

pdf(file.path(outdir, "Partial_RDA_Group_condition_InitialBW_baseplot.pdf"), width = 6, height = 5)
plot(rda_partial, scaling = 2, type = "n", main = "Partial RDA: Group | InitialBW")
points(
  rda_partial,
  display = "sites",
  scaling = 2,
  col = ifelse(meta_sub$Group == "Control", '#B22222', '#4169E1'),
  pch = 19,
  cex = 1.3
)
text(rda_partial, display = "cn", scaling = 2, col = "blue", cex = 0.9)
legend(
  "topright",
  legend = levels(meta_sub$Group),
  col = c('#B22222', '#4169E1'),
  pch = 19,
  bty = "n"
)
dev.off()

############################################################
# 14. DCA
############################################################

# DCA는 CCA/RDA 선택용 진단에 가까움
dca_res <- decorana(taxa_rel_filt)

cat("\n===== DCA result =====\n")
print(dca_res)

cat("\n===== DCA eigenvalues / decorana values =====\n")
print(eigenvals(dca_res))

pdf(file.path(outdir, "DCA_genus_baseplot.pdf"), width = 6, height = 5)
plot(dca_res, main = "DCA of genus profiles")
dev.off()

############################################################
# 15. CCA
############################################################

# 주의:
# 이 table은 raw count가 아니라 relative abundance이므로,
# 논문용 주분석은 CCA보다 Hellinger-RDA 또는 Bray-Curtis PCoA/dbRDA가 더 적절함.
# 아래 CCA는 비교/실습용으로 포함.

cca_res <- cca(
  taxa_rel_filt ~ InitialBW + Group,
  data = meta_sub
)

cat("\n===== CCA result: taxa_rel_filt ~ InitialBW + Group =====\n")
print(cca_res)

cat("\n===== CCA adjusted R2 =====\n")
print(RsquareAdj(cca_res))

cat("\n===== CCA ANOVA: full model =====\n")
print(anova(cca_res, permutations = 999))

cat("\n===== CCA ANOVA: by term =====\n")
print(anova(cca_res, by = "term", permutations = 999))

cat("\n===== CCA ANOVA: by margin =====\n")
print(anova(cca_res, by = "margin", permutations = 999))

cat("\n===== CCA ANOVA: by axis =====\n")
print(anova(cca_res, by = "axis", permutations = 999))

pdf(file.path(outdir, "CCA_Group_InitialBW_baseplot.pdf"), width = 6, height = 5)
plot(cca_res, scaling = 2, type = "n", main = "CCA: Group + InitialBW")
points(
  cca_res,
  display = "sites",
  scaling = 2,
  col = ifelse(meta_sub$Group == "Control", '#B22222', '#4169E1'),
  pch = 19,
  cex = 1.3
)
text(cca_res, display = "cn", scaling = 2, col = "blue", cex = 0.9)
legend(
  "topright",
  legend = levels(meta_sub$Group),
  col = c('#B22222', '#4169E1'),
  pch = 19,
  bty = "n"
)
dev.off()

############################################################
# 16. Bray-Curtis PCoA
############################################################

bc_dist <- vegdist(taxa_rel_filt, method = "bray")

pcoa_res <- cmdscale(
  bc_dist,
  eig = TRUE,
  k = 2
)

eig <- pcoa_res$eig
positive_eig <- eig[eig > 0]
pcoa_var <- positive_eig / sum(positive_eig) * 100

pcoa_df <- data.frame(
  SampleID = rownames(taxa_rel_filt),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  Group = meta_sub$Group,
  SamplingWeek = meta_sub$SamplingWeek,
  InitialBW = meta_sub$InitialBW
)

p_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Control" = '#B22222', "VNAM" = '#4169E1')) +
  labs(
    title = "Bray-Curtis PCoA of genus profiles",
    x = paste0("PCoA1 (", round(pcoa_var[1], 1), "%)"),
    y = paste0("PCoA2 (", round(pcoa_var[2], 1), "%)")
  ) +
  theme_classic(base_size = 14)

if (min(table(pcoa_df$Group)) >= 3) {
  p_pcoa <- p_pcoa + stat_ellipse(type = "t", linewidth = 0.7)
}

print(p_pcoa)

ggsave(
  filename = file.path(outdir, "PCoA_BrayCurtis_genus.pdf"),
  plot = p_pcoa,
  width = 5,
  height = 4
)

write.csv(
  pcoa_df,
  file = file.path(outdir, "PCoA_BrayCurtis_scores.csv"),
  row.names = FALSE
)

############################################################
# 17. PERMANOVA
############################################################

# 순차검정: InitialBW를 먼저 보정한 뒤 Group 설명력 확인
adonis_seq <- adonis2(
  bc_dist ~ InitialBW + Group,
  data = meta_sub,
  permutations = 999
)

cat("\n===== PERMANOVA: sequential test, InitialBW + Group =====\n")
print(adonis_seq)

# marginal test
adonis_margin <- adonis2(
  bc_dist ~ InitialBW + Group,
  data = meta_sub,
  permutations = 999,
  by = "margin"
)

cat("\n===== PERMANOVA: marginal test =====\n")
print(adonis_margin)

capture.output(
  adonis_seq,
  file = file.path(outdir, "PERMANOVA_sequential_InitialBW_Group.txt")
)

capture.output(
  adonis_margin,
  file = file.path(outdir, "PERMANOVA_margin_InitialBW_Group.txt")
)

############################################################
# 18. dbRDA / capscale
############################################################

# Bray-Curtis distance 기반 constrained ordination
dbrda_res <- capscale(
  taxa_rel_filt ~ InitialBW + Group,
  data = meta_sub,
  distance = "bray",
  add = TRUE
)

cat("\n===== dbRDA / capscale result =====\n")
print(dbrda_res)

cat("\n===== dbRDA adjusted R2 =====\n")
print(RsquareAdj(dbrda_res))

cat("\n===== dbRDA ANOVA: full model =====\n")
print(anova(dbrda_res, permutations = 999))

cat("\n===== dbRDA ANOVA: by term =====\n")
print(anova(dbrda_res, by = "term", permutations = 999))

cat("\n===== dbRDA ANOVA: by margin =====\n")
print(anova(dbrda_res, by = "margin", permutations = 999))

pdf(file.path(outdir, "dbRDA_BrayCurtis_Group_InitialBW_baseplot.pdf"), width = 6, height = 5)
plot(dbrda_res, scaling = 2, type = "n", main = "dbRDA: Bray-Curtis, Group + InitialBW")
points(
  dbrda_res,
  display = "sites",
  scaling = 2,
  col = ifelse(meta_sub$Group == "Control", '#B22222', '#4169E1'),
  pch = 19,
  cex = 1.3
)
text(dbrda_res, display = "cn", scaling = 2, col = "blue", cex = 0.9)
legend(
  "topright",
  legend = levels(meta_sub$Group),
  col = c('#B22222', '#4169E1'),
  pch = 19,
  bty = "n"
)
dev.off()

############################################################
# 19. Save processed matrices and metadata
############################################################

write.csv(
  taxa_rel,
  file = file.path(outdir, "taxa_relative_normalized_all.csv")
)

write.csv(
  taxa_rel_filt,
  file = file.path(outdir, "taxa_relative_filtered.csv")
)

write.csv(
  taxa_hell,
  file = file.path(outdir, "taxa_hellinger_filtered.csv")
)

write.csv(
  meta_sub,
  file = file.path(outdir, "metadata_matched.csv"),
  row.names = FALSE
)

############################################################
# 20. End
############################################################

cat("\n===== All analyses completed =====\n")
cat("Results saved in: ", outdir, "\n")