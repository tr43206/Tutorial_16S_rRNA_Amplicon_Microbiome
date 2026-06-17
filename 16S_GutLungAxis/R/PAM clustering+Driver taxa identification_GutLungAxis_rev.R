rm(list = ls())
setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")  #Edit

suppressPackageStartupMessages({
  library(cluster)      # pam, silhouette
  library(clusterSim)   # index.G1 (CH index)
  library(ade4)         # dudi.pca, dudi.pco, bca
  library(ggplot2)
  library(ggpubr)       # ggscatter (optional)
  library(factoextra)   # fviz_eig (optional)
})

#install.packages('ggfortify')
library("tidyverse")
library("qiime2R")
library("ggpubr")
library("ggfortify")
library("grid")
library("gridExtra")

#####################################################################################
### 원래 그룹대로 PCoA 그리기

metadata<-read_q2metadata("RUN45_Metadata_rev_GutLungAxis_phyloseq.txt")  #Edit
uwunifrac<-read_qza("core-metrics-results-20731/unweighted_unifrac_pcoa_results.qza")  #Edit
wunifrac<-read_qza("core-metrics-results-20731/weighted_unifrac_pcoa_results.qza")  #Edit


UnW<-uwunifrac$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata)

unw<-ggplot(UnW, aes(x=PC1, y=PC2, color=`Group`)) +  #Edit: color로 PCoA에 나타낼 그룹명 지정
  geom_point(alpha=0.8, size =2)+theme_q2r()+
  geom_text(data = subset(UnW, Group == "VNAM"),  #Edit: subset에 annotation 원하는 그룹명 지정
            aes(label = Group),
            vjust = -0.7,
            size = 2.5,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
unw


W<-wunifrac$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata)

w<-ggplot(W, aes(x=PC1, y=PC2, color=`Group`)) +  #Edit
  geom_point(alpha=0.8, size =2)+theme_q2r()+
  geom_text(data = subset(W, Group == "VNAM"),  #Edit
            aes(label = Group),
            vjust = -0.7,
            size = 2.5,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
w

#####################################################################################

infile <- 'rel_gg2_genus_table.tsv'

## Filtering params
MIN_PREV  <- 0.20  #Edit: 최소 출현 비율 (ex. 전체 샘플의 20% 이상에서)
MIN_REL   <- 1e-4  #Edit: 출현으로 간주할 최소 상대풍부도
MIN_TOTAL <- 1e-3  #Edit: 전체 합 기준 MIN_TOTAL 미만인 택사 제거

#####################################################################################
### Load data (taxa x samples; 1열=taxon ID)

## A1(1열 1행)만 공백으로 수정 후 읽기
L <- readLines(infile)
L[1] <- sub("^[^\t]+", "", L[1])  ## A1만 공백으로
writeLines(L, infile)

data <- read.table(infile, header=TRUE, row.names=1, sep="\t",
                   dec=".", check.names=FALSE, quote="", comment.char="")
data=data[-1,]  # 첫 행 전체 날림

# 행이름 중복 방지  #Edit: 여기부터 Noise filtering 단계까지는 선택사항
rownames(data) <- make.unique(rownames(data))
# 수치형 열만 유지 (문자열/주석열 방지)
data <- data[, vapply(data, is.numeric, logical(1)), drop=FALSE]

## 열 합=1 정규화
cs <- colSums(data)
if (any(cs == 0)) {
  keep_cols <- which(cs > 0)
  warning(sprintf("합계=0인 샘플 %d개 제거", sum(cs==0)))
  data <- data[, keep_cols, drop=FALSE]
  cs <- colSums(data)
}
data <- sweep(data, 2, cs, "/")

#####################################################################################
### Noise filtering (prevalence + total) -> 재정규화

filter_prevalence <- function(df, min_prev = 0.2, min_rel = 1e-4) {
  keep <- rowMeans(df > min_rel) >= min_prev
  df[keep, , drop = FALSE]
}
filter_total <- function(df, min_total = 1e-3) {
  keep <- rowSums(df) >= min_total
  df[keep, , drop = FALSE]
}

data <- filter_prevalence(data, MIN_PREV, MIN_REL)
data <- filter_total(data, MIN_TOTAL)

## 필터링 후 열 합=1 재정규화
data <- sweep(data, 2, colSums(data), "/")

#####################################################################################
### PAM clustering

pam.clustering <- function(d, k){
  cluster::pam(d, k, diss=TRUE)$clustering
}

pam_res <- pam(UnW[, c("PC1", "PC2")], k = 2)  #Edit: 원하는 클러스터 개수(k) 직접 지정

# 클러스터 결과를 UnW에 추가
UnW$Cluster <- as.factor(pam_res$clustering)  ## Cluster라는 컬럼에 1, 2로 저장됨 (factor로)
W$Cluster <- as.factor(pam_res$clustering)

cl_best <- pam_res$clustering

data2 <- data %>% t()
data3 <- data2[UnW$SampleID, ]
data4 <- data3 %>% t()

#####################################################################################
### Cluster <-> Sample mapping (k_best)

cluster_assignment <- data.frame(Subject = colnames(data4),
                                 Cluster = cl_best)
cluster_assignment <- cluster_assignment[order(cluster_assignment$Cluster), ]
print(cluster_assignment)

#####################################################################################

unw<-ggplot(UnW, aes(x=PC1, y=PC2, color=Cluster)) +  #Edit
  geom_point(alpha=0.8, size =2)+theme_q2r()+
  geom_text(data = subset(UnW, Group == "VNAM"),  #Edit
            aes(label = Group),
            vjust = -0.7,
            size = 2.5,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
unw


w<-ggplot(W, aes(x=PC1, y=PC2, color=`Cluster`)) +  #Edit
  geom_point(alpha=0.8, size =2)+theme_q2r()+
  geom_text(data = subset(W, Group == "VNAM"),  #Edit
            aes(label = Group),
            vjust = -0.7,
            size = 2.5,
            show.legend = FALSE) +
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
w

#####################################################################################
### Driver taxa by taxonomy level (Cluster 1 vs 2)

## cl_best <-> 샘플 정렬 일치시키기
if (is.null(names(cl_best))) {
  names(cl_best) <- colnames(data4)  ## 이름 없으면 순서대로 부여
} else {
  cl_best <- cl_best[colnames(data4)]  ## 이름 있으면 컬럼 순서에 맞춰 재정렬
}

## 계통 라벨 파서 (행이름 ex. d__Bacteria;p__Firmicutes;...;g__Bacteroides;...)
parse_lineage <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]  ## ; 기준으로 나누라는 뜻
  get <- function(prefix) {
    hit <- grep(paste0("^", prefix, "__"), parts)
    if (length(hit)) sub(paste0("^", prefix, "__"), "", parts[hit[1]]) else NA_character_
  }
  c(kingdom = get("d"), phylum = get("p"), class = get("c"),
    order   = get("o"), family = get("f"), genus = get("g"), species = get("s"))
}
tax_df <- as.data.frame(t(vapply(rownames(data4), parse_lineage,
                                 FUN.VALUE = rep(NA_character_, 7))),
                        stringsAsFactors = FALSE, row.names = rownames(data4))

## 한 레벨에서 드라이버 상위 n개 계산
top_drivers <- function(level, topn = 10, pseudo = 1e-8) {
  key <- tax_df[[level]]
  key[is.na(key) | key == ""] <- paste0("Unassigned_", level)
  # 같은 레벨 라벨끼리 합산(행 집계)
  agg <- rowsum(data4, group = key, reorder = FALSE)
  # 클러스터별 평균
  idx1 <- which(cl_best == 1)
  idx2 <- which(cl_best == 2)
  if (!length(idx1) || !length(idx2)) {
    return(data.frame())  ## 한쪽이 비면 빈 결과
  }
  m1 <- rowMeans(agg[, idx1, drop = FALSE])
  m2 <- rowMeans(agg[, idx2, drop = FALSE])
  diff   <- m2 - m1
  log2fc <- log2((m2 + pseudo) / (m1 + pseudo))
  out <- data.frame(level = rownames(agg),
                    mean_cl1 = m1, mean_cl2 = m2,
                    diff = diff, log2fc = log2fc,
                    row.names = NULL, check.names = FALSE)
  out <- out[order(abs(out$diff), decreasing = TRUE), ]  #Edit: diff or log2fc
  head(out, topn)
}

## 레벨별 표 생성 + 콘솔 출력 + CSV 저장
levels_vec <- c("genus","family","order","class","phylum")
drivers_list <- setNames(vector("list", length(levels_vec)), levels_vec)

out_dir <- getwd()  ## 현재 폴더에 저장
for (lev in levels_vec) {
  tab <- top_drivers(lev, topn = 10)
  drivers_list[[lev]] <- tab
  cat("\n=== Top drivers by", toupper(lev), "(|Δ mean| 상위 10) ===\n")
  if (nrow(tab)) print(tab, row.names = FALSE, digits = 4) else cat("No features.\n")
#  fn <- file.path(out_dir, sprintf("drivers_by_%s_top10.csv", lev))
#  utils::write.csv(tab, fn, row.names = FALSE)
#  cat(">> saved:", fn, "\n")
}  #Edit: 저장하려면 위의 fn부터 cat까지 주석 해제

# ## (옵션) 하나의 Excel로 저장하고 싶으면 주석 해제 (writexl 필요)
# if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
# writexl::write_xlsx(
#   c(list(clusters = data.frame(Sample = colnames(data), Cluster = cl_best)),
#     lapply(drivers_list, function(x) { if (is.null(x)) data.frame() else x })),
#   path = file.path(out_dir, "cluster_drivers_by_taxonomy.xlsx")
# )

## 각 csv(drivers_by_<level>_top10.csv)는 해당 레벨에서 클러스터2 - 클러스터1 평균 차이가 큰 순으로 정렬됨.
#  diff > 0 ⇒ 클러스터 2에서 높은 분류군
#  diff < 0 ⇒ 클러스터 1에서 높은 분류군
#  log2fc는 방향/크기 참고(0 분모 방지용 pseudo 포함)

#####################################################################################

## 4그룹(Control/VNAM × Cluster 1/2) 컬럼 추가
library(dplyr)

UnW <- UnW %>%
  mutate(
    Group4 = paste0(as.character(Group), "-", as.character(Cluster)),
    Group4 = factor(Group4, levels = c("Control-1","Control-2","VNAM-1","VNAM-2"))
  )

W <- W %>%
  mutate(
    Group4 = paste0(as.character(Group), "-", as.character(Cluster)),
    Group4 = factor(Group4, levels = c("Control-1","Control-2","VNAM-1","VNAM-2"))
  )

#####################################################################################

unw<-ggplot(UnW, aes(x=PC1, y=PC2, color=Group4)) +  #Edit
  geom_point(alpha=0.8, size =2)+theme_q2r()+
#  geom_text(data = subset(UnW, Group == "VNAM"),  #Edit
#            aes(label = Group),
#            vjust = -0.7,
#            size = 2.5,
#            show.legend = FALSE) +  #Edit: 지금 주석처리 되어있는 5줄 주석처리 풀면 원하는 그룹 annotation 할 수 있음
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
unw


w<-ggplot(W, aes(x=PC1, y=PC2, color=Group4)) +  #Edit
  geom_point(alpha=0.8, size =2)+theme_q2r()+
#  geom_text(data = subset(W, Group == "VNAM"),  #Edit
#            aes(label = Group),
#            vjust = -0.7,
#            size = 2.5,
#            show.legend = FALSE) +
  scale_colour_manual(values = c("#5975A4", "#CC8963", "#5F9E6E","#B55D60"))+
  xlab("") +
  ylab("") +
  font("xlab", size = 15, face="bold")+
  font("ylab", size = 15, face="bold")+
  font("xy.text", size = 10, face="bold")+
  geom_vline(xintercept=0.0, colour="lightgrey", size=0.7) +
  geom_hline(yintercept=0.0, colour="lightgrey",  size=0.7)+
  theme(legend.position = "right")+
  stat_ellipse()
w

#####################################################################################

## Control-1 / VNAM-1만 필터링
UnW1 <- UnW %>% filter(Group4 %in% c("Control-1","VNAM-1"))
W1   <- W   %>% filter(Group4 %in% c("Control-1","VNAM-1"))

## 색 2개만 지정 (Control-1, VNAM-1 순서 고정)
cols1 <- c("Control-1" = "#5975A4", "VNAM-1" = "#5F9E6E")

unw1 <- ggplot(UnW1, aes(x = PC1, y = PC2, color = Group4)) +
  geom_point(alpha = 0.8, size = 2) + theme_q2r() +
  scale_colour_manual(values = cols1, drop = FALSE) +
  xlab("") + ylab("") +
  font("xlab", size = 15, face = "bold") +
  font("ylab", size = 15, face = "bold") +
  font("xy.text", size = 10, face = "bold") +
  geom_vline(xintercept = 0.0, colour = "lightgrey", size = 0.7) +
  geom_hline(yintercept = 0.0, colour = "lightgrey", size = 0.7) +
  theme(legend.position = "right") +
  stat_ellipse()

unw1


w1 <- ggplot(W1, aes(x = PC1, y = PC2, color = Group4)) +
  geom_point(alpha = 0.8, size = 2) + theme_q2r() +
  scale_colour_manual(values = cols1, drop = FALSE) +
  xlab("") + ylab("") +
  font("xlab", size = 15, face = "bold") +
  font("ylab", size = 15, face = "bold") +
  font("xy.text", size = 10, face = "bold") +
  geom_vline(xintercept = 0.0, colour = "lightgrey", size = 0.7) +
  geom_hline(yintercept = 0.0, colour = "lightgrey", size = 0.7) +
  theme(legend.position = "right") +
  stat_ellipse()

w1

#####################################################################################

## Control-2 / VNAM-2만 필터링
UnW2 <- UnW %>% filter(Group4 %in% c("Control-2","VNAM-2"))
W2   <- W   %>% filter(Group4 %in% c("Control-2","VNAM-2"))

## 색 2개만 지정 (Control-2, VNAM-2 순서 고정)
cols2 <- c("Control-2" = "#CC8963", "VNAM-2" = "#B55D60")

unw2 <- ggplot(UnW2, aes(x = PC1, y = PC2, color = Group4)) +
  geom_point(alpha = 0.8, size = 2) + theme_q2r() +
  scale_colour_manual(values = cols2, drop = FALSE) +
  xlab("") + ylab("") +
  font("xlab", size = 15, face = "bold") +
  font("ylab", size = 15, face = "bold") +
  font("xy.text", size = 10, face = "bold") +
  geom_vline(xintercept = 0.0, colour = "lightgrey", size = 0.7) +
  geom_hline(yintercept = 0.0, colour = "lightgrey", size = 0.7) +
  theme(legend.position = "right") +
  stat_ellipse()

unw2


w2 <- ggplot(W2, aes(x = PC1, y = PC2, color = Group4)) +
  geom_point(alpha = 0.8, size = 2) + theme_q2r() +
  scale_colour_manual(values = cols2, drop = FALSE) +
  xlab("") + ylab("") +
  font("xlab", size = 15, face = "bold") +
  font("ylab", size = 15, face = "bold") +
  font("xy.text", size = 10, face = "bold") +
  geom_vline(xintercept = 0.0, colour = "lightgrey", size = 0.7) +
  geom_hline(yintercept = 0.0, colour = "lightgrey", size = 0.7) +
  theme(legend.position = "right") +
  stat_ellipse()

w2

#####################################################################################
#####################################################################################
### distance 계산

## qiime2 view에서 받은 distance.tsv 읽기 (행/열 = SampleID)
dist_df <- read.delim("beta_div/uw_distance.tsv",
                      header = TRUE,
                      sep = '\t',
                      check.names = FALSE) %>%
  dplyr::select(SubjectID1, SubjectID2, Distance)


## UnW에서 SampleID - Group4 매핑만 추출 (중복 있으면 제거)
map <- UnW %>%
  dplyr::select(SampleID, Group4) %>%
  distinct()

dist_g <- dist_df %>%
  left_join(map, by = c("SubjectID1" = "SampleID")) %>% rename(Group4_1 = Group4) %>%
  left_join(map, by = c("SubjectID2" = "SampleID")) %>% rename(Group4_2 = Group4) %>%
  filter(!is.na(Group4_1), !is.na(Group4_2)) %>%
  mutate(
    pair_type = if_else(Group4_1 == Group4_2, "within", "between"),
    g_min = pmin(as.character(Group4_1), as.character(Group4_2)),
    g_max = pmax(as.character(Group4_1), as.character(Group4_2)),
    pair_group = if_else(pair_type == "within",
                         paste0(g_min, "_within"),
                         paste0(g_min, "_vs_", g_max))
  )

dist_summary <- dist_g %>%
  group_by(pair_group, pair_type) %>%
  summarise(n = n(),
            mean_dist = mean(Distance, na.rm=TRUE),
            sd_dist   = sd(Distance, na.rm=TRUE),
            .groups="drop") %>%
  arrange(pair_type, pair_group)

dist_summary

#####################################################################################
### PERMANOVA - 그룹 간 유의한 쌍이 하나 이상 있는지 확인

library(vegan)

## dist_df(long) -> dist 객체 만들기
dl <- dist_df %>% distinct(SubjectID1, SubjectID2, .keep_all=TRUE)
samp <- sort(unique(c(dl$SubjectID1, dl$SubjectID2)))

mat <- matrix(0, length(samp), length(samp), dimnames=list(samp, samp))
i <- match(dl$SubjectID1, samp)
j <- match(dl$SubjectID2, samp)
mat[cbind(i,j)] <- dl$Distance
mat[cbind(j,i)] <- dl$Distance
diag(mat) <- 0
d <- as.dist(mat)

## Group4를 UnW에서 바로 매칭 (파일 불필요)
grp4 <- setNames(as.character(UnW$Group4), UnW$SampleID)[samp]
keep <- !is.na(grp4)

d2 <- as.dist(mat[keep, keep])
meta_use <- data.frame(Group4 = factor(grp4[keep]))

set.seed(1)
adonis2(d2 ~ Group4, data = meta_use, permutations = 999)

#####################################################################################
### Pairwise p-val 계산

lv <- levels(meta_use$Group4)
pairs <- combn(lv, 2, simplify = FALSE)

Dmat <- as.matrix(d2)

out <- lapply(pairs, function(p){
  idx <- meta_use$Group4 %in% p
  d_sub <- as.dist(Dmat[idx, idx])
  dat_sub <- droplevels(meta_use[idx, , drop=FALSE])
  
  fit <- adonis2(d_sub ~ Group4, data = dat_sub, permutations = 999)
  data.frame(g1=p[1], g2=p[2],
             R2=fit$R2[1], F=fit$F[1], p=fit$`Pr(>F)`[1])
})

pairwise_tbl <- do.call(rbind, out)
pairwise_tbl$q <- p.adjust(pairwise_tbl$p, method="BH")
pairwise_tbl

#####################################################################################

