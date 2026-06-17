rm(list = ls())

#BiocManager::install("phyloseq")
#install.packages("ggpicrust2")
#BiocManager::install("ALDEx2")
library(phyloseq)   # 마이크로바이옴 데이터 분석 및 시각화
library(tidyverse)  # R의 데이터 핸들링 
library(readr)      # 파일 읽어오기 
library(ggpicrust2) # PICRUSt2결과 처리 및 시각화
library(ALDEx2)     # Aldex2 분석
library(mutoss)     # FDR(BKY) 분석
library(pheatmap)   # heatmap 시각화

#####################################################################################

setwd('C:/Users/tr432/Downloads/filter_GutLungAxis_Week0')

df <- read.csv("gg2_genus_table.tsv", sep = '\t', row.names = 1,
               skip = 1, check.names = FALSE, comment.char = "")
head(df)
dim(df)

#####################################################################################

meta0 <- read.csv("../filter_GutLungAxis_Control/maaslin3_RUN45_Metadata_GutLungAxis.txt", sep = '\t',
                      check.names = FALSE)
colnames(meta0)[1] <- "SampleID"
meta <- meta0[meta0$Group %in% c("Control", "VNAM"), ]
meta <- meta[meta$SampleID %in% colnames(df), ]  ## df에 실제로 존재하는 샘플만 다시 필터
df2 <- df[, meta$SampleID]
dim(df2)

#####################################################################################
##**## (예시)
# meta$Group <- factor(meta$Group, levels = c("Control", "VNAM"))
# meta$Sex   <- factor(meta$Sex)        # 범주형 공변량 예시
# meta$Age   <- as.numeric(meta$Age)    # 연속형 공변량 예시
##**##

Group <- factor(meta$Group)
#BM_wbc   <- as.numeric(meta$BM_wbc)
#b.w <- as.numeric(meta$b.w)
Initial_b.w <- as.numeric(meta$Initial_b.w)

Level <- levels(Group)


covariates <- data.frame(
  Group = meta$Group,
#  BM_wbc   = meta$BM_wbc#,
#  b.w   = meta$b.w,
  Initial_b.w   = meta$Initial_b.w
  # 필요한 공변량 더 추가 가능
)
rownames(covariates) <- rownames(meta)  # 샘플 순서 맞춰두기


design <- model.matrix(~ Group + Initial_b.w, data = covariates)
head(design)


ALDEx2_object <- ALDEx2::aldex.clr(
  reads       = round(df2),
  conds       = design,        # <- 여기 중요: model.matrix
  mc.samples  = 128,           # 원래 쓰던 값
  denom       = "all"          # 원래 설정
)


glm_res <- ALDEx2::aldex.glm(ALDEx2_object, verbose = TRUE)


# 1) GLM effect 전체 계산
glm_eff_list <- ALDEx2::aldex.glm.effect(ALDEx2_object)

# 어떤 contrast 들이 있는지 확인
names(glm_eff_list)
# 예: "(Intercept)", "GroupVNAM", "BM_wbc", ...

# 2) 관심 contrast 선택 (예: GroupVNAM)
eff_group <- glm_eff_list[["GroupVNAM"]]   # data.frame, rownames = feature
head(eff_group)
# 여기 column 중에 'effect' 가 ALDEx2 effect 와 같은 정의


colnames(glm_res)


# 1) 관심 있는 열 이름 자동으로 찾기
p_col   <- grep("^GroupVNAM:pval$",     colnames(glm_res), value = TRUE)
q_col   <- grep("^GroupVNAM:pval.padj$|^GroupVNAM:BH$", 
                colnames(glm_res), value = TRUE)

# 2) feature / p값 / q값 / effect만 모아서 새로운 데이터프레임 생성
glm_group <- data.frame(
  feature = rownames(glm_res),
  pval    = glm_res[[p_col]],
  qval    = glm_res[[q_col]],
  row.names = NULL
)

head(glm_group)


p_values_df <- data.frame(
  feature = glm_group$feature,
  group1  = Level[1],                 # 예: "Control"
  group2  = Level[2],                 # 예: "VNAM"
  p_values = glm_group$pval,
  stringsAsFactors = FALSE
)

dim(p_values_df)
head(p_values_df)


## (선택1) 패키지 사용해서 BKY method FDR

p_wi <- p_values_df %>% as.data.frame()
# p-value 벡터 확보
pv <- p_wi$p_values
# 결과 컬럼을 먼저 만들어서 "전체 feature 유지"
p_wi$q_bky <- NA_real_
# NA 아닌 것만 BKY 계산 후, 해당 위치에만 채우기
keep <- !is.na(pv)
bky_wi <- mutoss::two.stage(pValues = pv[keep], alpha = 0.05)
p_wi$q_bky[keep] <- bky_wi$adjPValues

bky_wi <- mutoss::two.stage(
  pValues = p_wi$p_values,
  alpha   = 0.05
)

p_wi$q_bky <- bky_wi$adjPValues

tax_AlDEx_result <- p_wi %>%
  dplyr::rename(p_adjust = q_bky)


##**##
## (선택2) 패키지 없이 BKY method FDR

bky_qvals <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  q <- rep(NA_real_, length(p))
  
  keep <- !is.na(p) & is.finite(p)
  pv <- p[keep]
  m  <- length(pv)
  if (m == 0) return(q)
  
  # Stage 1: BH at alpha/(1+alpha)
  alpha1 <- alpha / (1 + alpha)
  R1 <- sum(p.adjust(pv, method = "BH") <= alpha1)
  
  m0_hat <- max(1L, m - R1)
  # BKY q-values: scaled BH
  qv <- pmin(1, (m0_hat / m) * p.adjust(pv, method = "BH"))
  
  q[keep] <- qv
  q
}

p_wi <- p_values_df %>% as.data.frame()
p_wi$q_bky <- bky_qvals(p_wi$p_values, alpha = 0.05)

tax_AlDEx_result <- dplyr::rename(p_wi, p_adjust = q_bky)
##**##


eff_df <- eff_group %>%
  rownames_to_column("feature") %>%
  dplyr::select(feature, effect)

# 4) feature 기준으로 tax_AlDEx_result에 effect 열 추가
#    (기존에 effect 열이 있다면 먼저 지우고 덮어쓰기)
tax_AlDEx_result <- tax_AlDEx_result %>%
  dplyr::select(-any_of("effect")) %>%       # 있으면 삭제, 없으면 무시
  left_join(eff_df, by = "feature")


dim(tax_AlDEx_result)
head(tax_AlDEx_result)

#####################################################################################

# p-value값 보정이 유의한 것만 추출
#tax_AlDEx_result2 <- tax_AlDEx_result[tax_AlDEx_result$p_adjust < 0.05, ]
#dim(tax_AlDEx_result2)

#outdir <- 'C:/Users/tr432/Downloads/Figures_GutLungAxis/no0/Taxa/ALDEx2'
write_csv(tax_AlDEx_result, file.path("AlDEx2_Group+Initial_b.w_BKY.csv"))

#####################################################################################
