#Reference: "Microbiome confounders and quantitative profiling challenge predicted microbial targets in colorectal cancer development", Nature medicine, 2024.
# https://doi.org/10.1038/s41591-024-02963-2

################################################################################

rm(list = ls())

library(vegan)
library(tidyverse)

setwd('C:/Users/tr432/Downloads/filter_GutLungAxis_Week2')

################################################################################
# =========================
# 0. Input files
# =========================
feature_file <- 'gg2_genus_table.tsv'
metadata_file <- 'C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt'

out_dir <- 'covariate_dbRDA_results'
dir.create(out_dir, showWarnings = FALSE)

################################################################################
# =========================
# 1. Read genus table
# =========================
# 첫 줄: "# Constructed from biom file" 제거
genus_raw <- read.delim(
  feature_file,
  sep = "\t",
  skip = 1,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(genus_raw)[1] <- "Taxon"

# genus level만 유지할지 여부
# TRUE  = g__가 있는 assigned genus만 사용
# FALSE = unclassified row까지 포함
KEEP_ONLY_ASSIGNED_GENUS <- TRUE

if (KEEP_ONLY_ASSIGNED_GENUS) {
  genus_raw <- genus_raw %>%
    filter(grepl(";g__[^;]+$", Taxon))
}

# taxonomy 이름 짧게 정리
get_last_taxon <- function(x) {
  sapply(strsplit(x, ";"), function(z) {
    z <- z[z != ""]
    tail(z, 1)
  })
}

taxon_short <- get_last_taxon(genus_raw$Taxon)
taxon_short <- make.unique(taxon_short)

count_mat <- genus_raw %>%
  select(-Taxon) %>%
  as.data.frame()

rownames(count_mat) <- taxon_short

# taxa × sample -> sample × taxa
abund_count <- as.data.frame(t(count_mat))
abund_count <- abund_count %>%
  mutate(across(everything(), as.numeric))

################################################################################
# =========================
# 2. Read metadata
# =========================
meta <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA", "N/A", "NaN")
)

colnames(meta)[colnames(meta) == "#SampleID"] <- "SampleID"
rownames(meta) <- meta$SampleID

################################################################################
# =========================
# 3. Match samples
# =========================
common_samples <- intersect(rownames(abund_count), rownames(meta))

abund_count <- abund_count[common_samples, , drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

cat("Matched samples:", length(common_samples), "\n")
print(common_samples)

################################################################################
# =========================
# 4. Feature filtering
# =========================
MIN_PREVALENCE <- 2  #Edit
MIN_TOTAL_COUNT <- 10  #Edit

keep_features <- colSums(abund_count > 0, na.rm = TRUE) >= MIN_PREVALENCE &
  colSums(abund_count, na.rm = TRUE) >= MIN_TOTAL_COUNT

abund_count <- abund_count[, keep_features, drop = FALSE]

cat("Remaining taxa:", ncol(abund_count), "\n")

################################################################################
# =========================
# 5. Convert to relative abundance
# =========================
# Bray-Curtis dbRDA에서는 raw count보다 relative abundance 권장
abund_rel <- sweep(abund_count, 1, rowSums(abund_count, na.rm = TRUE), "/")
abund_rel[is.na(abund_rel)] <- 0

################################################################################
# =========================
# 6. Metadata covariates
# =========================
# Group은 main variable.
# covariate 후보는 별도로 둠.

#** 주의할 점
# 서로 거의 같은 batch/cohort 정보를 담고 있는 covariates를 동시에 넣으면 collinearity가 생김 (ex. SacrificeDate, ExpNo, SamplingDate).
# Collinearity (공선성): 회귀분석에서 예측 변수(독립 변수)들 사이에 강한 선형 관계가 있어, 서로의 변동을 중복 설명하는 현상.

main_variable <- 'Group'

covariate_candidates <- c(
  'ExpNo',
  'SamplingDate',
  'SacrificeDate',
  'Initial_b.w',
  'SequencingDepth'
)

# MGX_Reads는 shotgun read 수라서 16S genus table 분석에서는 기본 제외.
# 넣고 싶으면 아래 주석 해제.
# covariate_candidates <- c(covariate_candidates, "MGX_Reads")

screen_variables <- c(main_variable, covariate_candidates)
screen_variables <- screen_variables[screen_variables %in% colnames(meta)]

meta_cov <- meta[, screen_variables, drop = FALSE]

# numeric variables
numeric_vars <- c(
  'Initial_b.w',
  'SequencingDepth',
  'MGX_Reads',
  'b.w',
  'BM_wbc',
  'BM_whole',
  'Spleen',
  'Spleen/b.w'
)

for (v in intersect(numeric_vars, colnames(meta_cov))) {
  meta_cov[[v]] <- as.numeric(meta_cov[[v]])
}

# factor variables
factor_vars <- setdiff(colnames(meta_cov), numeric_vars)

for (v in factor_vars) {
  meta_cov[[v]] <- as.factor(meta_cov[[v]])
}

# 값이 1종류뿐인 변수 제거
keep_cov <- sapply(meta_cov, function(x) {
  length(unique(na.omit(x))) >= 2
})

meta_cov <- meta_cov[, keep_cov, drop = FALSE]

cat("Variables used for screening:\n")
print(colnames(meta_cov))

################################################################################
# =========================
# 7. Univariate dbRDA

# Nature 논문식 univariate dbRDA
# : 각 metadata variable이 Bray–Curtis community variation을 얼마나 설명하는지 평가함.
# =========================
run_univariate_dbrda <- function(abund_mat, metadata, permutations = 9999) {
  
  results <- list()
  
  for (var in colnames(metadata)) {
    
    tmp_meta <- metadata[, var, drop = FALSE]
    colnames(tmp_meta) <- "var"
    
    idx <- complete.cases(tmp_meta)
    
    x <- abund_mat[idx, , drop = FALSE]
    m <- tmp_meta[idx, , drop = FALSE]
    
    # factor level 정리
    if (is.factor(m$var)) {
      m$var <- droplevels(m$var)
    }
    
    # 값이 1종류면 skip
    if (length(unique(na.omit(m$var))) < 2) next
    
    # taxa filtering after NA removal
    x <- x[, colSums(x, na.rm = TRUE) > 0, drop = FALSE]
    
    if (nrow(x) < 4) next
    if (ncol(x) < 2) next
    
    fit <- tryCatch(
      capscale(
        x ~ var,
        data = m,
        distance = "bray",
        add = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    an <- tryCatch(
      anova.cca(fit, permutations = permutations),
      error = function(e) NULL
    )
    
    if (is.null(an)) next
    
    r2 <- RsquareAdj(fit)
    
    results[[var]] <- data.frame(
      Variable = var,
      N = nrow(x),
      R2 = r2$r.squared,
      Adj_R2 = r2$adj.r.squared,
      F_value = an$F[1],
      P_value = an$`Pr(>F)`[1]
    )
  }
  
  res <- bind_rows(results)
  
  res <- res %>%
    mutate(
      Q_value_BH = p.adjust(P_value, method = "BH")
    ) %>%
    arrange(Q_value_BH, desc(R2))
  
  return(res)
}

dbrda_res <- run_univariate_dbrda(
  abund_mat = abund_rel,
  metadata = meta_cov,
  permutations = 9999
)

write.csv(
  dbrda_res,
  file.path(out_dir, "01_univariate_dbRDA_covariates.csv"),
  row.names = FALSE
)

dbrda_res

################################################################################
# =========================
# 8. Select significant covariates

# Group 제외 후 non-redundant covariate 선정
# : Group은 treatment effect를 보는 main variable이므로, confounder 후보 선정에서는 제외함.
# =========================
ALPHA <- 0.05

sig_covariates <- dbrda_res %>%
  filter(
    Variable != main_variable,
    Q_value_BH < ALPHA
  ) %>%
  pull(Variable)

cat("Significant covariates after BH correction:\n")
print(sig_covariates)

################################################################################
# =========================
# 9. Stepwise dbRDA for non-redundant covariates
# =========================
run_stepwise_dbrda <- function(abund_mat, metadata, variables, permutations = 9999) {
  
  if (length(variables) == 0) {
    message("No significant covariates. Stepwise dbRDA skipped.")
    return(NULL)
  }
  
  # =========================
  # 1. 변수 제한
  # =========================
  variables <- variables[variables %in% colnames(metadata)]
  
  if (length(variables) == 0) {
    message("No available covariates in metadata.")
    return(NULL)
  }
  
  m <- metadata[, variables, drop = FALSE]
  idx <- complete.cases(m)
  
  x <- abund_mat[idx, , drop = FALSE]
  m <- m[idx, , drop = FALSE]
  m <- droplevels(m)
  
  # all-zero taxa 제거
  x <- x[, colSums(x, na.rm = TRUE) > 0, drop = FALSE]
  
  # 값이 1종류뿐인 변수 제거
  keep_vars <- sapply(m, function(v) {
    length(unique(na.omit(v))) >= 2
  })
  
  m <- m[, keep_vars, drop = FALSE]
  
  if (ncol(m) == 0) {
    message("No valid covariates after filtering constant variables.")
    return(NULL)
  }
  
  # =========================
  # 2. 변수명 안전하게 변경
  # =========================
  name_map <- data.frame(
    original = colnames(m),
    safe = make.names(colnames(m), unique = TRUE),
    stringsAsFactors = FALSE
  )
  
  colnames(m) <- name_map$safe
  
  # =========================
  # 3. ordiR2step 평가 문제 방지
  # =========================
  assign(".dbrda_x", x, envir = .GlobalEnv)
  assign(".dbrda_m", m, envir = .GlobalEnv)
  
  on.exit({
    rm(".dbrda_x", envir = .GlobalEnv)
    rm(".dbrda_m", envir = .GlobalEnv)
  }, add = TRUE)
  
  # =========================
  # 4. dbRDA null / full model
  # =========================
  mod0 <- capscale(
    .dbrda_x ~ 1,
    data = .dbrda_m,
    distance = "bray",
    add = TRUE
  )
  
  mod_full <- capscale(
    .dbrda_x ~ .,
    data = .dbrda_m,
    distance = "bray",
    add = TRUE
  )
  
  # =========================
  # 5. Forward stepwise dbRDA
  # =========================
  step_mod <- ordiR2step(
    mod0,
    scope = formula(mod_full),
    direction = "forward",
    permutations = permutations,
    Pin = 0.05,
    R2scope = TRUE,
    trace = TRUE
  )
  
  selected_safe <- attr(terms(step_mod), "term.labels")
  
  selected_covariates <- name_map %>%
    filter(safe %in% selected_safe)
  
  # =========================
  # 6. 결과 정리
  # =========================
  anova_overall <- anova.cca(step_mod, permutations = permutations)
  anova_terms <- anova.cca(step_mod, by = "term", permutations = permutations)
  r2_res <- RsquareAdj(step_mod)
  
  return(list(
    model = step_mod,
    name_map = name_map,
    selected_covariates = selected_covariates,
    anova_overall = anova_overall,
    anova_terms = anova_terms,
    r2 = r2_res
  ))
}

step_res <- run_stepwise_dbrda(
  abund_mat = abund_rel,
  metadata = meta_cov,
  variables = sig_covariates,
  permutations = 9999
)

if (!is.null(step_res)) {
  
  write.csv(
    step_res$selected_covariates,
    file.path(out_dir, "02_stepwise_dbRDA_selected_covariates.csv"),
    row.names = FALSE
  )
  
  print(step_res$selected_covariates)
  print(step_res$r2)
  print(step_res$anova_terms)
}

################################################################################
# =========================
# 10. Adjusted PERMANOVA

# Group effect를 selected covariates로 보정해서 PERMANOVA
# : Group 효과가 covariate 보정 후에도 community-level에서 남는지 봄.
# =========================
selected_covariates <- character(0)

if (!is.null(step_res)) {
  selected_covariates <- step_res$selected_covariates$original
}

adonis_vars <- c(main_variable, selected_covariates)
adonis_vars <- adonis_vars[adonis_vars %in% colnames(meta_cov)]

m_adonis <- meta_cov[, adonis_vars, drop = FALSE]
idx <- complete.cases(m_adonis)

x_adonis <- abund_rel[idx, , drop = FALSE]
m_adonis <- droplevels(m_adonis[idx, , drop = FALSE])

bray <- vegdist(x_adonis, method = "bray")

formula_adonis <- as.formula(
  paste("bray ~", paste(adonis_vars, collapse = " + "))
)

adonis_res <- adonis2(
  formula_adonis,
  data = m_adonis,
  permutations = 9999,
  by = "margin"
)

write.csv(
  as.data.frame(adonis_res),
  file.path(out_dir, "03_adjusted_PERMANOVA_Group_effect.csv")
)

adonis_res

################################################################################
# =========================
# 11. Unadjusted taxon-level test

# Taxon-level: 보정 전 Kruskal–Wallis
# : 각 genus가 Control group vs. Treatment group에서 다른지 먼저 봄.
# =========================
group_vec <- meta[rownames(abund_rel), main_variable]
group_vec <- as.factor(group_vec)

kw_res <- lapply(colnames(abund_rel), function(taxon) {
  
  df <- data.frame(
    Group = group_vec,
    Abundance = abund_rel[, taxon]
  )
  
  df <- df[complete.cases(df), ]
  
  if (length(unique(df$Group)) < 2) return(NULL)
  if (length(unique(df$Abundance)) < 2) return(NULL)
  
  test <- kruskal.test(Abundance ~ Group, data = df)
  
  data.frame(
    Taxon = taxon,
    P_value = test$p.value,
    Statistic = as.numeric(test$statistic)
  )
})

kw_res <- bind_rows(kw_res) %>%
  mutate(Q_value_BH = p.adjust(P_value, method = "BH")) %>%
  arrange(Q_value_BH)

write.csv(
  kw_res,
  file.path(out_dir, "04_taxa_unadjusted_Kruskal.csv"),
  row.names = FALSE
)

kw_res

################################################################################
# =========================
# 12. Covariate-adjusted taxon-level test

# Taxon-level: selected covariates 보정 후 Group 효과
# : Rank-transformed abundance ~ Group + covariates 형태로 봄.
#   샘플 수가 적으면 logistic GLM보다 이 방식이 덜 불안정함.
# =========================
run_adjusted_lm_nested <- function(taxon, covariates) {
  
  df <- data.frame(
    Group = meta[rownames(abund_rel), main_variable],
    Abundance = abund_rel[, taxon],
    meta[rownames(abund_rel), covariates, drop = FALSE]
  )
  
  df <- df[complete.cases(df), ]
  
  if (nrow(df) < 6) return(NULL)
  if (length(unique(df$Group)) < 2) return(NULL)
  if (length(unique(df$Abundance)) < 2) return(NULL)
  
  df$Group <- as.factor(df$Group)
  df$Abundance_rank <- rank(df$Abundance, ties.method = "average")
  
  # covariate 처리
  cov_terms <- c()
  
  for (cv in covariates) {
    
    if (is.numeric(df[[cv]])) {
      new_cv <- paste0(cv, "_rank")
      df[[new_cv]] <- rank(df[[cv]], ties.method = "average")
      cov_terms <- c(cov_terms, paste0("`", new_cv, "`"))
      
    } else {
      df[[cv]] <- as.factor(df[[cv]])
      
      # level이 1개면 제외
      if (length(unique(df[[cv]])) >= 2) {
        cov_terms <- c(cov_terms, paste0("`", cv, "`"))
      }
    }
  }
  
  # reduced model: covariates only
  if (length(cov_terms) == 0) {
    reduced_formula <- as.formula("Abundance_rank ~ 1")
  } else {
    reduced_formula <- as.formula(
      paste("Abundance_rank ~", paste(cov_terms, collapse = " + "))
    )
  }
  
  # full model: covariates + Group
  full_formula <- as.formula(
    paste(
      "Abundance_rank ~",
      paste(c(cov_terms, "Group"), collapse = " + ")
    )
  )
  
  fit_reduced <- tryCatch(lm(reduced_formula, data = df), error = function(e) NULL)
  fit_full <- tryCatch(lm(full_formula, data = df), error = function(e) NULL)
  
  if (is.null(fit_reduced) | is.null(fit_full)) return(NULL)
  
  test <- anova(fit_reduced, fit_full)
  
  data.frame(
    Taxon = taxon,
    N = nrow(df),
    P_value_Group_adjusted = test$`Pr(>F)`[2],
    F_value_Group = test$F[2]
  )
}


selected_covariates <- 'SacrificeDate'  #Edit

adj_taxa_res <- lapply(colnames(abund_rel), function(taxon) {
  run_adjusted_lm_nested(taxon, selected_covariates)
}) %>%
  bind_rows()

adj_taxa_res <- adj_taxa_res %>%
  mutate(Q_value_BH_adjusted = p.adjust(P_value_Group_adjusted, method = "BH")) %>%
  arrange(Q_value_BH_adjusted)

write.csv(
  adj_taxa_res,
  file.path(out_dir, "05_taxa_covariate_adjusted_lm.csv"),
  row.names = FALSE
)

adj_taxa_res

################################################################################
# =========================
# 13. Before vs after adjustment

# 보정 전/후 비교표
# =========================
compare_res <- kw_res %>%
  select(
    Taxon,
    KW_P = P_value,
    KW_Q = Q_value_BH
  ) %>%
  left_join(
    adj_taxa_res %>%
      select(
        Taxon,
        Adjusted_P = P_value_Group_adjusted,
        Adjusted_Q = Q_value_BH_adjusted
      ),
    by = "Taxon"
  ) %>%
  mutate(
    Before_adjustment = KW_Q < 0.05,
    After_adjustment = Adjusted_Q < 0.05,
    Interpretation = case_when(
      Before_adjustment & After_adjustment ~ "Robust after covariate adjustment",
      Before_adjustment & !After_adjustment ~ "Potentially confounded or weakened after adjustment",
      !Before_adjustment & After_adjustment ~ "Only significant after adjustment",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(Adjusted_Q)

write.csv(
  compare_res,
  file.path(out_dir, "06_taxa_before_after_covariate_adjustment.csv"),
  row.names = FALSE
)

compare_res

################################################################################

## 결과 해석 기준

# 1. 전체 community variation에 영향을 주는 변수
dbrda_res
# -> Q_value_BH < 0.05인 covariates 선택

# 2. Group 제외 후 non-redundant covariate
step_res$selected_covariates
# -> stepwise dbRDA 단계에서 covariate-1을 이미 보정한 상태에서 covariate-2가 추가 설명력(R2)을 가지는가?를 판단함 (ex. covariate-1 + 2 vs. covariate-1 설명력 비교). 만약 추가 설명력(R2)을 가지는 2번째 변수가 없다면 변수들 중에서 설명력(R2)이 가장 높은 변수를 결과로 출력.

# 3. covariate 보정 후에도 Treatment group effect가 남는지
adonis_res
# -> "Group" vs. "Group + covariate-1" 비교. Covariate-1을 보정해도 treatment group의 microbiome composition 차이가 유의하게 남는지를 판단.

# 4. genus별 보정 전/후 유의성 변화
head(compare_res, 5)
sum(compare_res$Before_adjustment == TRUE, na.rm = TRUE)  ## 보정 전 유의한 taxa 개수
sum(compare_res$After_adjustment == TRUE, na.rm = TRUE)  ## 보정 후 유의한 taxa 개수
## 보정 후 유의해지지 않은 taxa
compare_res %>%
  filter(Before_adjustment == TRUE & After_adjustment == FALSE)
## 보정 후 유의해진 taxa
compare_res %>%
  filter(Before_adjustment == FALSE & After_adjustment == TRUE)
## 특정 taxa만 확인
compare_res %>%
  filter(grepl("Phocaeicola", Taxon))  #Edit
# -> KW_Q vs. Adjucted_Q

################################################################################
################################################################################

library(tidyverse)
library(ggplot2)

out_dir <- "covariate_dbRDA_figures"
dir.create(out_dir, showWarnings = FALSE)

################################################################################
# =========================
# 1) Univariate dbRDA result

# Covariate effect size plot
# =========================
plot_dbrda <- dbrda_res %>%
  mutate(
    Variable = factor(Variable, levels = rev(Variable)),
    R2_percent = R2 * 100,
    Adj_R2_percent = Adj_R2 * 100,
    Significant = ifelse(Q_value_BH < 0.05, "FDR < 0.05", "n.s.")
  )

p1 <- ggplot(plot_dbrda, aes(x = R2_percent, y = Variable, fill = Significant)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(round(R2_percent, 1), "%")),
    hjust = -0.15,
    size = 4
  ) +
  scale_fill_manual(values = c("FDR < 0.05" = "firebrick", "n.s." = "gray80")) +
  labs(
    x = "Explained variation by univariate dbRDA (%)",
    y = NULL,
    fill = NULL,
    title = "Individual covariate effects on genus-level microbiome variation"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  coord_cartesian(xlim = c(0, max(plot_dbrda$R2_percent, na.rm = TRUE) * 1.25))

p1

ggsave(
  filename = file.path(out_dir, "Fig_covariate_univariate_dbRDA.png"),
  plot = p1,
  width = 8,
  height = 6,
  dpi = 500
)

################################################################################
# =========================
# 2) Adjusted PERMANOVA plot

# PERMANOVA 결과 barplot
# -> Group effect가 보정 후에도 남는지 표시
# =========================
adonis_df <- as.data.frame(adonis_res) %>%
  rownames_to_column("Variable") %>%
  filter(!Variable %in% c("Residual", "Total")) %>%
  mutate(
    R2_percent = R2 * 100,
    P_label = case_when(
      `Pr(>F)` < 0.001 ~ "***",
      `Pr(>F)` < 0.01  ~ "**",
      `Pr(>F)` < 0.05  ~ "*",
      TRUE ~ "n.s."
    ),
    Variable = factor(Variable, levels = rev(Variable))
  )

p2 <- ggplot(adonis_df, aes(x = R2_percent, y = Variable)) +
  geom_col(width = 0.75, fill = "gray70", color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = paste0(round(R2_percent, 1), "% ", P_label)),
    hjust = -0.15,
    size = 4.5
  ) +
  labs(
    x = "Marginal R² from adjusted PERMANOVA (%)",
    y = NULL,
    title = "VNAM group effect after covariate adjustment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  coord_cartesian(xlim = c(0, max(adonis_df$R2_percent, na.rm = TRUE) * 1.25))

p2

ggsave(
  filename = file.path(out_dir, "Fig_adjusted_PERMANOVA_Group_effect.png"),
  plot = p2,
  width = 6.5,
  height = 3,
  dpi = 500
)

################################################################################
# =========================
# 3) Summary count plot

# Robust / Lost / Gained taxa 개수 plot
# =========================
summary_count <- compare_res %>%
  mutate(
    Category = case_when(
      Before_adjustment & After_adjustment ~ "Robust",
      Before_adjustment & !After_adjustment ~ "Lost after adjustment",
      !Before_adjustment & After_adjustment ~ "Gained after adjustment",
      TRUE ~ "Not significant"
    )
  ) %>%
  count(Category) %>%
  mutate(
    Category = factor(
      Category,
      levels = c("Robust", "Lost after adjustment", "Gained after adjustment", "Not significant")
    )
  )

p3 <- ggplot(summary_count, aes(x = Category, y = n, fill = Category)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n), vjust = -0.4, size = 5) +
  scale_fill_manual(
    values = c(
      "Robust" = "firebrick",
      "Lost after adjustment" = "gray50",
      "Gained after adjustment" = "royalblue",
      "Not significant" = "gray85"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of taxa",
    title = "Effect of covariate adjustment on genus-level associations"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(0, max(summary_count$n) * 1.1))

p3

ggsave(
  filename = file.path(out_dir, "Fig_taxa_adjustment_summary_count.png"),
  plot = p3,
  width = 6.5,
  height = 6.5,
  dpi = 500
)

################################################################################
# =========================
# 4) Before/After q-value comparison

# 보정 전/후 q-value lollipop plot
# -> 각 taxon의 보정 전후 유의성 강도를 보여줌.
# =========================
TOP_N <- 40

plot_lollipop <- compare_res %>%
  mutate(
    before_score = -log10(KW_Q),
    after_score = -log10(Adjusted_Q),
    max_score = pmax(before_score, after_score, na.rm = TRUE),
    Status = case_when(
      Before_adjustment & After_adjustment ~ "Robust",
      Before_adjustment & !After_adjustment ~ "Lost after adjustment",
      !Before_adjustment & After_adjustment ~ "Gained after adjustment",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(desc(max_score)) %>%
  slice_head(n = TOP_N) %>%
  mutate(Taxon = factor(Taxon, levels = rev(Taxon)))

plot_lollipop_long <- plot_lollipop %>%
  select(Taxon, before_score, after_score, Status) %>%
  pivot_longer(
    cols = c(before_score, after_score),
    names_to = "Analysis",
    values_to = "minus_log10_q"
  ) %>%
  mutate(
    Analysis = recode(
      Analysis,
      "before_score" = "Before adjustment",
      "after_score" = "After adjustment"
    ),
    Analysis = factor(Analysis, levels = c("Before adjustment", "After adjustment"))
  )

p4 <- ggplot(plot_lollipop, aes(y = Taxon)) +
  geom_segment(
    aes(x = before_score, xend = after_score, yend = Taxon),
    color = "gray60",
    linewidth = 0.6
  ) +
  geom_point(
    data = plot_lollipop_long,
    aes(x = minus_log10_q, color = Analysis),
    size = 2.8
  ) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Before adjustment" = "gray30",
      "After adjustment" = "firebrick"
    )
  ) +
  labs(
    x = expression(-log[10]("FDR q-value")),
    y = NULL,
    color = NULL,
    title = "Before vs after covariate adjustment"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9.5, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

p4

ggsave(
  filename = file.path(out_dir, "Fig_taxa_before_after_qvalue_lollipop.png"),
  plot = p4,
  width = 7,
  height = 9,
  dpi = 500
)

################################################################################
# =========================
# 5) Phocaeicola_A before/after adjustment
# =========================
phocaeicola_res <- compare_res %>%
  filter(grepl("Phocaeicola", Taxon)) %>%
  select(Taxon, KW_Q, Adjusted_Q) %>%
  pivot_longer(
    cols = c(KW_Q, Adjusted_Q),
    names_to = "Analysis",
    values_to = "Q_value"
  ) %>%
  mutate(
    Analysis = recode(
      Analysis,
      "KW_Q" = "Before adjustment",
      "Adjusted_Q" = "After adjustment"
    ),
    minus_log10_q = -log10(Q_value),
    Significant = ifelse(Q_value < 0.05, "FDR < 0.05", "n.s.")
  )

p5 <- ggplot(phocaeicola_res, aes(x = Analysis, y = minus_log10_q, fill = Significant)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5) +
  geom_text(
    aes(label = paste0("q = ", signif(Q_value, 3))),
    vjust = -0.4,
    size = 4.5
  ) +
  scale_fill_manual(values = c("FDR < 0.05" = "firebrick", "n.s." = "gray80")) +
  labs(
    x = NULL,
    y = expression(-log[10]("FDR q-value")),
    fill = NULL,
    title = "Phocaeicola_A association before and after adjustment"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  coord_cartesian(
    ylim = c(0, max(phocaeicola_res$minus_log10_q, na.rm = TRUE) * 1.15)
  )

p5

ggsave(
  filename = file.path(out_dir, "Fig_Phocaeicola_before_after_adjustment.png"),
  plot = p5,
  width = 6,
  height = 6,
  dpi = 500
)

################################################################################
