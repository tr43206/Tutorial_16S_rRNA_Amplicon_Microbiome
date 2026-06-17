rm(list = ls())

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biobakery/maaslin3")
library(maaslin3)
library(mutoss)
library(dplyr)
library(ggplot2)

#####################################################################################
### Microbiome association detection with MaAsLin 3

setwd("C:/Users/tr432/Desktop/Roh/_Study_Design/25.05.02~_Gut-Lung_Axis/Exp1+1-S/MGX/04.humann.out/pathabundance_relab_split")

### 상대적 풍부도 데이터
taxa_table <- read.csv("pathabundance_relab_stratified_maaslin3.txt",
                       sep = '\t',
                       row.names = 1,
#                       skip = 1,
                       check.names = FALSE, comment.char = "")
taxa_table <- t(taxa_table)


metadata <- read.csv("../metadata.csv", row.names = 1)


##**## (예시)
# meta$Group <- factor(meta$Group, levels = c("Control", "VNAM"))
# meta$Sex   <- factor(meta$Sex)        # 범주형 공변량 예시
# meta$Age   <- as.numeric(meta$Age)    # 연속형 공변량 예시
##**##

metadata$Group <- factor(metadata$Group, levels = c('Control', 'VNAM'))
#metadata$Initial_b.w   <- as.numeric(metadata$Initial_b.w)

colnames(taxa_table) <- gsub("\\.", "-", colnames(taxa_table))

taxa_table[1:5, 1:5]
#metadata[1:5, 1:5]
metadata

#####################################################################################

set.seed(1)
fit_out <- maaslin3(
  input_data    = taxa_table,
  input_metadata = metadata,
  output        = "maaslin3_Group_BH_stratified",  #Edit
  formula       = "~ Group",  #Edit
  normalization = "TSS",
  transform     = "LOG",  ## "LOG"는 base 2 log / "PLOG"는 pseudo-log
#  zero_threshold = 0,      ## 26.01.02 추가 - transform = "PLOG" 기준
#  evaluate_only  = "abundance",   ## 26.01.02 추가 - transform = "PLOG" 기준
#  warn_prevalence = FALSE,  ## 26.01.02 추가 - transform = "PLOG" 기준
  correction    = "BH",
  augment       = TRUE,
  standardize   = TRUE,
  max_significance = 0.1,
  median_comparison_abundance  = TRUE,  ## 상대적 풍부도엔 써도되지만 절대 풍부도엔 쓰면 안됨
  median_comparison_prevalence = TRUE,
  cores        = 1
)

## * 샘플 수가 적은데 (약 20개 내외) CageNo 종류가 많아서 오류가 났을 것으로 예상_25.12.15

#####################################################################################

## Abundance - BKY FDR
res_abund <- fit_out$fit_data_abundance$results %>%
  dplyr::select(-qval_individual, -qval_joint)

## individual p-values에 대한 BKY
keep_abund_ind <- !is.na(res_abund$pval_individual)
bky_abund_ind  <- mutoss::two.stage(res_abund$pval_individual[keep_abund_ind], alpha = 0.05)
res_abund$q_abund_ind_bky <- NA_real_
res_abund$q_abund_ind_bky[keep_abund_ind] <- bky_abund_ind$adjPValues


##**##
## (선택) 만약 유의한게 없어서 오류가 난다면 아래 부분으로 대체 (BKY 보정 직접 보정)
pv <- suppressWarnings(as.numeric(res_abund$pval_individual))
keep <- is.finite(pv)

bky_qvalues <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  m <- length(p)
  o <- order(p)
  ro <- order(o)
  ps <- p[o]
  
  # stage 1: BH at alpha' = alpha/(1+alpha)
  alpha1 <- alpha / (1 + alpha)
  crit1  <- (1:m) * alpha1 / m
  idx1   <- which(ps <= crit1)
  r1     <- if (length(idx1) == 0) 0 else max(idx1)
  
  # m0-hat = m - r1  (BKY two-step estimate)
  m0 <- m - r1
  if (m0 < 1) m0 <- 1  # 극단 케이스(거의 없음) 안전장치
  
  # stage 2: BH-style adjusted p-values with m0 in denominator
  qs <- (m0 * ps) / (1:m)
  qs <- rev(cummin(rev(qs)))
  qs[qs > 1] <- 1
  
  qs[ro]
}

res_abund$q_abund_ind_bky <- NA_real_
if (any(keep)) {
  # 1) mutoss 시도
  out <- try(mutoss::two.stage(pv[keep], alpha = 0.05), silent = TRUE)
  
  adj <- NULL
  if (!inherits(out, "try-error")) {
    adj <- out[["adjPValues"]]
    if (!is.null(adj)) adj <- as.numeric(adj)
  }
  
  # 2) mutoss가 이상하면 직접 계산으로 폴백
  if (is.null(adj) || length(adj) != sum(keep)) {
    adj <- bky_qvalues(pv[keep], alpha = 0.05)
  }
  
  res_abund$q_abund_ind_bky[keep] <- adj
}
##**##

write.csv(res_abund, "maaslin3_Group_abund_BKY_stratified.csv", row.names = FALSE)  #Edit

#####################################################################################

## Prevalence - BKY FDR
res_prev <- fit_out$fit_data_prevalence$results %>%
  dplyr::select(-qval_individual, -qval_joint)

## individual p-values에 대한 BKY
keep_prev_ind <- !is.na(res_prev$pval_individual)
bky_prev_ind  <- mutoss::two.stage(res_prev$pval_individual[keep_prev_ind], alpha = 0.05)
res_prev$q_prev_ind_bky <- NA_real_
res_prev$q_prev_ind_bky[keep_prev_ind] <- bky_prev_ind$adjPValues


##**##
## (선택) 만약 유의한게 없어서 오류가 난다면 아래 부분으로 대체 (BKY 보정 직접 보정)
pv <- suppressWarnings(as.numeric(res_prev$pval_individual))
keep <- is.finite(pv)

bky_qvalues <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  m <- length(p)
  o <- order(p)
  ro <- order(o)
  ps <- p[o]

  # stage 1: BH at alpha' = alpha/(1+alpha)
  alpha1 <- alpha / (1 + alpha)
  crit1  <- (1:m) * alpha1 / m
  idx1   <- which(ps <= crit1)
  r1     <- if (length(idx1) == 0) 0 else max(idx1)

  # m0-hat = m - r1  (BKY two-step estimate)
  m0 <- m - r1
  if (m0 < 1) m0 <- 1  # 극단 케이스(거의 없음) 안전장치

  # stage 2: BH-style adjusted p-values with m0 in denominator
  qs <- (m0 * ps) / (1:m)
  qs <- rev(cummin(rev(qs)))
  qs[qs > 1] <- 1

  qs[ro]
}

res_prev$q_prev_ind_bky <- NA_real_
if (any(keep)) {
# 1) mutoss 시도
  out <- try(mutoss::two.stage(pv[keep], alpha = 0.05), silent = TRUE)

  adj <- NULL
  if (!inherits(out, "try-error")) {
    adj <- out[["adjPValues"]]
    if (!is.null(adj)) adj <- as.numeric(adj)
  }

# 2) mutoss가 이상하면 직접 계산으로 폴백
  if (is.null(adj) || length(adj) != sum(keep)) {
    adj <- bky_qvalues(pv[keep], alpha = 0.05)
  }

  res_prev$q_prev_ind_bky[keep] <- adj
}
##**##

write.csv(res_prev, "maaslin3_Group_prev_BKY_stratified.csv", row.names = FALSE)  #Edit

#####################################################################################
