# 0) RStudio 재시작 직후 실행

# 1) lock 폴더 제거
#unlink("C:/Users/tr432/AppData/Local/R/win-library/4.5/00LOCK",
#       recursive = TRUE, force = TRUE)

# 2) 기존 패키지 제거
#remove.packages("ANCOMBC")
#remove.packages("CVXR")

# 3) remotes 설치
#install.packages("remotes")

# 4) CVXR를 1.0-15로 설치
#remotes::install_version(
#  "CVXR",
#  version = "1.8.1",
#  repos = "https://cran.rstudio.com/"
#)

# 5) ANCOMBC 재설치
#BiocManager::install("ANCOMBC", force = TRUE, ask = FALSE, update = FALSE)

# 6) 확인
#packageVersion("CVXR")
#packageVersion("ANCOMBC")

#library(ANCOMBC)

## =========================
## 0. Install / Load
## =========================

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")

rm(list = ls())

library(ANCOMBC)
library(dplyr)
library(mutoss)

setwd('C:/Users/tr432/Downloads/filter_GutLungAxis_Week0')

#####################################################################################
## =========================
## 1. Input files
## =========================

## 1) 읽기
otu0  <- read.delim("gg2_genus_table.tsv", sep="\t", header=TRUE, skip = 1,
                    row.names=1, check.names=FALSE, comment.char="", quote="")
meta0 <- read.delim("C:/Users/tr432/Downloads/filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt", sep="\t", header=TRUE,
                    row.names=1, check.names=FALSE, comment.char="", quote="")

otu0  <- as.data.frame(otu0)
meta0 <- as.data.frame(meta0)

## 2) ID 정리 함수 (BOM/공백/제어문자 제거)
clean_id <- function(x){
  x <- as.character(x)
  x <- sub("^\ufeff", "", x)          # BOM 제거
  x <- gsub("[[:cntrl:]]", "", x)     # 숨은 제어문자 제거
  x <- trimws(x)                      # 앞뒤 공백 제거
  x
}

## 3) metadata에서 "샘플ID 컬럼" 자동 탐지 (가장 많이 매칭되는 컬럼 선택)
otu_col_ids <- clean_id(colnames(otu0))
otu_row_ids <- clean_id(rownames(otu0))

score_col <- sapply(names(meta0), function(nm){
  v <- meta0[[nm]]
  if (!(is.character(v) || is.factor(v))) return(0)
  length(intersect(clean_id(v), otu_col_ids))
})
score_row <- sapply(names(meta0), function(nm){
  v <- meta0[[nm]]
  if (!(is.character(v) || is.factor(v))) return(0)
  length(intersect(clean_id(v), otu_row_ids))
})

best_col <- names(which.max(score_col))
best_row <- names(which.max(score_row))

## 현재 meta rownames도 후보로 점수 계산
score_rn_col <- length(intersect(clean_id(rownames(meta0)), otu_col_ids))
score_rn_row <- length(intersect(clean_id(rownames(meta0)), otu_row_ids))

## 4) 어떤 방향/어떤 컬럼이 가장 잘 맞는지 결정
cand <- c(
  rn_col = score_rn_col,
  rn_row = score_rn_row,
  bestcol_col = max(score_col),
  bestrow_row = max(score_row)
)
print(sort(cand, decreasing=TRUE))

## 5) metadata rownames 확정 (샘플ID 컬럼이 더 잘 맞으면 그걸 rownames로 설정)
if (max(score_col) > score_rn_col || max(score_row) > score_rn_row) {
  # 둘 중 더 큰 쪽 선택
  if (max(score_col) >= max(score_row)) {
    message("metadata 샘플ID 컬럼으로 추정: ", best_col, " (otu colnames와 매칭)")
    rownames(meta0) <- clean_id(meta0[[best_col]])
    meta0[[best_col]] <- NULL
  } else {
    message("metadata 샘플ID 컬럼으로 추정: ", best_row, " (otu rownames와 매칭)")
    rownames(meta0) <- clean_id(meta0[[best_row]])
    meta0[[best_row]] <- NULL
  }
} else {
  # meta rownames가 이미 샘플ID로 보임
  rownames(meta0) <- clean_id(rownames(meta0))
}

## 6) otu 방향 결정 (샘플이 col인지 row인지)
# 정리된 meta 샘플ID
sid <- rownames(meta0)

col_ok <- length(intersect(sid, otu_col_ids))
row_ok <- length(intersect(sid, otu_row_ids))

message("매칭 개수: meta vs otu colnames = ", col_ok, " / meta vs otu rownames = ", row_ok)

if (col_ok >= row_ok) {
  otu <- as.matrix(otu0)
  colnames(otu) <- otu_col_ids
  rownames(otu) <- clean_id(rownames(otu0))
} else {
  otu <- t(as.matrix(otu0))
  colnames(otu) <- otu_row_ids
  rownames(otu) <- clean_id(colnames(otu0))
}

## 7) 교집합으로 정렬 (완전매칭 강제 X)
sid <- clean_id(rownames(meta0))
rownames(meta0) <- sid

common <- intersect(sid, colnames(otu))
message("최종 공통 샘플 수 = ", length(common))

if (length(common) < 3) stop("공통 샘플이 너무 적습니다. 샘플ID 포맷이 크게 다릅니다(파일 확인 필요).")

otu  <- otu[, common, drop=FALSE]
meta <- meta0[common, , drop=FALSE]

## 8) 중복 샘플ID 체크
if (anyDuplicated(colnames(otu))) stop("otu 샘플ID 중복 존재")
if (anyDuplicated(rownames(meta))) stop("metadata 샘플ID 중복 존재")

## 결과 확인
head(colnames(otu))
head(rownames(meta))

#####################################################################################
## 3.2 Run ancombc2 function

meta$Group <- factor(meta$Group)
meta$`Initial_b.w` <- as.numeric(meta$`Initial_b.w`)
meta$`SequencingDepth` <- as.numeric(meta$`SequencingDepth`)
meta$`MGX_Reads` <- as.numeric(meta$`MGX_Reads`)


## =========================
## 6. Run ANCOM-BC2
## =========================
set.seed(123)

out <- ancombc2(
  data = otu,
  meta_data = meta,
  fix_formula = "Group + Initial_b.w",
#  rand_formula = NULL,
  group = "Group",
  prv_cut = 0,
  lib_cut = 0,
#  s0_perc = 0.05,
  p_adj_method = "BH",
  pseudo = 1,
  pseudo_sens = TRUE,  #Edit
  struc_zero = TRUE,   #Edit
  neg_lb = TRUE,       #Edit
  alpha = 0.05,
#  n_cl = 2,
  verbose = TRUE,
  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE
)

res_prim <- out$res
#res_pair <- out$res_pair

#out$zero_ind

######################################################################################################################

## Structural Zero 비율/개수 확인 (struc_zero = TRUE 일 때만 작동함)
#* -> 만약 비율이 20-30퍼 이상이라면 struc_zero = TRUE 로 켜는걸 권장.
#* -> 그래도 유의한 feature 수가 많지 않다면 neg_lb = TRUE 도 켜는걸 권장 (negative lower bound).
## --> Setting neg_lb = TRUE indicates that you are using both criteria stated in section 3.2 of ANCOM-II to detect structural zeros; otherwise, the algorithm will only use the equation 1 in section 3.2 for declaring structural zeros. Generally, it is recommended to set neg_lb = TRUE when the sample size per group is relatively large (e.g. > 30).

zi <- out$zero_ind
dim(zi)
head(zi)

# taxon별: 어떤 그룹에서라도 structural zero로 뜨면 TRUE
tax_any <- apply(zi, 1, function(x) any(x == TRUE, na.rm = TRUE))
mean(tax_any)*100          # 전체 taxon 중 structural zero taxon 비율
sum(tax_any)               # 개수

######################################################################################################################

## =========================
## 7. BKY (two-stage) adjustment on p-value columns
## =========================

bky_adjust <- function(p, alpha=0.05) {
  p <- suppressWarnings(as.numeric(p))
  outq <- rep(NA_real_, length(p))
  
  ok <- is.finite(p) & !is.na(p) & p >= 0 & p <= 1
  m  <- sum(ok)
  
  if (m == 0) return(outq)
  if (m == 1) { outq[ok] <- p[ok]; return(outq) }  # 유효 p 1개면 그대로 둠(원하면 1로 바꿔도 됨)
  
  tmp <- try(mutoss::two.stage(p[ok], alpha = alpha), silent = TRUE)
  
  if (inherits(tmp, "try-error") || is.null(tmp$adjp) || length(tmp$adjp) != m) {
    # BKY 실패 시 안전한 fallback(BH)
    outq[ok] <- p.adjust(p[ok], method = "BH")
  } else {
    outq[ok] <- tmp$adjp
  }
  outq
}

add_bky <- function(df, alpha=0.05) {
  if (is.null(df) || NROW(df) == 0) return(df)
  
  # grep() 말고 grepl() 사용 (논리 벡터)
  p_cols <- names(df)[
    grepl("^p_", names(df)) |
      grepl("pval|p_val|pvalue", names(df), ignore.case = TRUE)
  ]
  
  for (pc in p_cols) {
    df[[paste0("q_BKY__", pc)]] <- bky_adjust(df[[pc]], alpha = alpha)
  }
  df
}

res_prim_bky <- add_bky(out$res, alpha=0.05)
#res_pair_bky <- if (!is.null(out$res_pair)) add_bky(out$res_pair, alpha=0.05) else NULL


## =========================
## 9. Save / View
## =========================
write.csv(res_prim_bky, "ANCOMBC2_Group+Initial_b.w+SequencingDepth_primary_BKY_pseudo1.csv", row.names=FALSE)
#write.csv(res_pair_bky, "ANCOMBC2_pairwise_BKY.csv", row.names=FALSE)

#####################################################################################
