rm(list = ls())

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# =========================================================
# 0) 경로/파일
# =========================================================
setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Control")
meta_path <- "maaslin3_RUN45_Metadata_GutLungAxis.txt"

# =========================================================
# 1) 메타데이터 로드 (문자 그대로 읽기)
# =========================================================
meta_raw <- read.delim(
  meta_path,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# =========================================================
# 2) 기본 정리/형변환 (필요한 컬럼만 안전하게)
# =========================================================
meta <- meta_raw %>%
  mutate(
    # Group
    Group = factor(Group, levels = c("Control", "VNAM")),
    trt   = ifelse(Group == "VNAM", 1, 0),
    
    # Age: "8W" -> 8
    Age_wk = suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(Age)))),
    
    # Sex / CageNo
    Sex    = factor(Sex),
    CageNo = factor(CageNo),
    
    # 숫자형: N/A -> NA
    Initial_b.w = suppressWarnings(as.numeric(Initial_b.w)),
    `b.w`       = suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(`b.w`)))),
    SequencingDepth = suppressWarnings(as.numeric(SequencingDepth))
  )

# =========================================================
# 3) 분석 단위 선택 (예: Week2만)
# =========================================================
WK <- "Week0"
meta0 <- meta %>%
  filter(SamplingWeek == WK) %>%
  filter(!is.na(Group))

cat("\n===== Group counts (", WK, ") =====\n", sep = "")
print(table(meta0$Group, useNA = "ifany"))

# =========================================================
# 4) 공변량 분포 점검: (1) SMD/비율차 테이블 + (2) 플롯
#    - 매칭/PSM 없이 "불균형(=교란 가능성)"만 확인
# =========================================================
covars <- c("ExpNo", "SamplingDate", "Initial_b.w", "SacrificeDate", "SequencingDepth")

# 연속형 SMD
smd_cont <- function(x, g){
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) return(NA_real_)
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  if (all(is.na(x1)) || all(is.na(x2))) return(NA_real_)
  m1 <- mean(x1, na.rm=TRUE); m2 <- mean(x2, na.rm=TRUE)
  s1 <- sd(x1, na.rm=TRUE);   s2 <- sd(x2, na.rm=TRUE)
  sp <- sqrt((s1^2 + s2^2)/2)
  if (is.na(sp) || sp == 0) return(NA_real_)
  (m2 - m1) / sp
}

# 범주형: 그룹별 분포 차이의 최대값(간단 지표)
propdiff_cat <- function(x, g){
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) return(NA_real_)
  tb <- prop.table(table(x, g), margin=2)
  if (nrow(tb) == 0) return(NA_real_)
  max(abs(tb[,2] - tb[,1]), na.rm=TRUE)
}

# ---- Binary SMD (0/1) ----
smd_binary <- function(z, g){
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) return(NA_real_)
  z1 <- z[g == levels(g)[1]]
  z2 <- z[g == levels(g)[2]]
  p1 <- mean(z1, na.rm = TRUE)
  p2 <- mean(z2, na.rm = TRUE)
  sp <- sqrt((p1*(1-p1) + p2*(1-p2))/2)
  if (is.na(sp) || sp == 0) return(NA_real_)
  (p2 - p1) / sp
}

# ---- Categorical SMD: 레벨별 더미 SMD 계산 후 요약 ----
smd_cat_summary <- function(x, g){
  x <- droplevels(factor(x))
  levs <- levels(x)
  if (length(levs) < 2) {
    return(list(value = NA_real_, details = NULL))
  }
  smds <- sapply(levs, function(L){
    z <- as.integer(x == L)
    smd_binary(z, g)
  })
  details <- data.frame(level = levs, smd = smds, abs_smd = abs(smds)) %>%
    arrange(desc(abs_smd))
  list(
    value   = details$abs_smd[1],   # 요약: max(|SMD|)
    details = details
  )
}


# 지표 계산: numeric=연속형 SMD, factor=범주형 max(|SMD|) 1개로 요약
cat_details_list <- list()

# 지표 계산
bal_tbl <- lapply(covars, function(v){
  x <- meta0[[v]]
  
  if (is.numeric(x)) {
    data.frame(var=v, type="numeric", metric="SMD", value=smd_cont(x, meta0$Group))
  } else {
    out <- smd_cat_summary(x, meta0$Group)
    # 레벨별 SMD 상세 저장(원하면 나중에 확인/저장 가능)
    if (!is.null(out$details)) cat_details_list[[v]] <<- out$details
    data.frame(var=v, type="factor",  metric="SMD", value=out$value)
  }
}) %>% bind_rows() %>%
  mutate(abs_value = abs(value)) %>%
  arrange(desc(abs_value))

cat("\n===== Imbalance summary (no matching) =====\n")
print(bal_tbl)

# (선택) 범주형 레벨별 SMD 상세를 파일로 저장하고 싶으면
# for (nm in names(cat_details_list)) {
#   write.table(cat_details_list[[nm]],
#               file = paste0("SMD_levels_", WK, "_", nm, ".tsv"),
#               sep="\t", quote=FALSE, row.names=FALSE)
# }

######################################################################################

##**##
# covars 안에서 실제 존재하는 컬럼만
covars <- intersect(covars, names(meta0))

# 숫자형/비숫자형 분리
num_covars <- covars[sapply(meta0[covars], is.numeric)]
cat_covars <- covars[!sapply(meta0[covars], is.numeric)]

# (1) 숫자형용 long
meta_long_num <- meta0 %>%
  select(Group, any_of(num_covars)) %>%
  pivot_longer(-Group, names_to="var", values_to="val")

# (2) 범주형(문자/팩터)용 long
meta_long_cat <- meta0 %>%
  mutate(across(any_of(cat_covars), as.factor)) %>%
  select(Group, any_of(cat_covars)) %>%
  pivot_longer(-Group, names_to="var", values_to="val")
##**##

######################################################################################
## Numerical

##**##
# 그룹 개수
n_group <- nlevels(meta_long_num$Group)

label_df_num <- bal_tbl %>%  ## bal_tbl: var, metric, value 포함
  filter(var %in% num_covars) %>%
  mutate(label = paste0(metric, " = ", round(value, 3)),
         x_pos = (n_group + 1) / 2)
##**##


# (A) 연속형 분포 플롯
p_num <- ggplot(meta_long_num, aes(x=Group, y=val, fill=Group)) +
  geom_violin(alpha=0.25, trim=FALSE) +
  geom_boxplot(width=0.15, outlier.shape=NA, alpha=0.6) +
  geom_jitter(width=0.1, size=1, alpha=0.8, color="black") +
  facet_wrap(~var, scales="free_y") +
#  scale_y_continuous(limits = c(19,23),
#                     breaks = seq(19,23, by = 0.5),
#                     expand = c(0,0)) +
  theme_classic(base_size = 14) +
  theme(legend.position="none",
        axis.title = element_text(size=12, face="bold"),
        axis.text  = element_text(size=10),
        strip.text = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, face = "bold", hjust=0)#,
#        panel.border = element_rect(colour = "black",
#                                    fill = NA,
#                                    linewidth = 1.5)
        ) +
  labs(title="Covariate distributions - numerical factors", x=NULL, y=NULL) +
  geom_text(data = label_df_num,
            aes(x = x_pos, y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 1.5,size = 4)

print(p_num)

ggsave("C:/Users/tr432/Downloads/WeeksMerged_16sSeq/Covariates_adequacy/covariates_numerical_w0.png",
       p_num, width = 7, height = 5, dpi = 500)

######################################################################################
## Categorical

##**##
# 변수별 레벨 수 계산
level_info <- meta_long_cat %>%
  group_by(var) %>%
  summarise(n_level = n_distinct(val), .groups="drop")

label_df_cat <- bal_tbl %>%
  filter(var %in% cat_covars) %>%
  left_join(level_info, by="var") %>%
  mutate(label = paste0(metric, " = ", round(value, 3)),
         x_pos = (n_level + 1) / 2)
##**##


# (B) 범주형 비율 플롯
p_cat <- ggplot(meta_long_cat, aes(x=val, fill=Group)) +
  geom_bar(position="fill") +
  facet_wrap(~var, scales="free_x") +
  scale_y_continuous(limits = c(0,1.08),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size=12, face="bold"),
        axis.text  = element_text(size=10),
        strip.text = element_text(size=12, face="bold"),
        plot.title = element_text(size=12, face="bold", hjust=0)#,
#        panel.border = element_rect(colour = "black",
#                                    fill = NA,
#                                    linewidth = 1)
        ) +
  labs(title="Covariate proportions - categorical factors", x=NULL, y="Proportion") +
  geom_text(data = label_df_cat,
            aes(x = x_pos, y = 1.04, label = label),
            inherit.aes = FALSE,
            size = 4)

print(p_cat)

ggsave("C:/Users/tr432/Downloads/WeeksMerged_16sSeq/Covariates_adequacy/covariates_categorical_찐.png",
       p_cat, width = 8, height = 5, dpi = 500)
