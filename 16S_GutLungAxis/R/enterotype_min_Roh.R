### Packages

suppressPackageStartupMessages({
  library(cluster)      # pam, silhouette
  library(clusterSim)   # index.G1 (CH index)
  library(ade4)         # dudi.pca, dudi.pco, bca
  library(ggplot2)
  library(ggpubr)       # ggscatter (optional)
  library(factoextra)   # fviz_eig (optional)
})

#####################################################################################

### Paths / Params

setwd('C:/Users/tr432/Downloads')
infile <- 'relative-genus_SML.tsv'

## Filtering params
MIN_PREV  <- 0.20    # 최소 출현 비율 (예: 전체 샘플의 20% 이상에서)
MIN_REL   <- 1e-4    # 출현으로 간주할 최소 상대풍부도
MIN_TOTAL <- 1e-3    # 전체 합 기준(0.1%) 미만인 택사 제거

## Clustering params
K_RANGE   <- 2:10
PSEUDO    <- 1e-6    # JSD용 pseudocount

#####################################################################################

### Load data (taxa x samples; 1열=taxon ID)

## A1(헤더 좌상단)만 공백으로 수정 후 읽기
L <- readLines(infile)
L[1] <- sub("^[^\t]+", "", L[1])   # 첫 필드(A1)만 공백으로
writeLines(L, infile)

data <- read.table(infile, header=TRUE, row.names=1, sep="\t",
                   dec=".", check.names=FALSE, quote="", comment.char="")
data=data[-1,]  # 첫 행 전체 날림

# 행이름 중복 방지
rownames(data) <- make.unique(rownames(data))
# 수치형 열만 유지 (문자열/주석열 방지)
data <- data[, vapply(data, is.numeric, logical(1)), drop=FALSE]

## 열 합=1 정규화 (JSD는 확률분포 가정)
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

### JSD distance (sqrt(JS) metric; log base 2)

dist.JSD <- function(mat, pseudocount=1e-6){
  mat <- as.matrix(mat)
  mat[mat == 0] <- pseudocount
  KLD <- function(x, y) sum(x * log2(x / y))
  JSD <- function(x, y) {
    m <- 0.5 * (x + y)
    0.5 * KLD(x, m) + 0.5 * KLD(y, m)
  }
  n <- ncol(mat); cn <- colnames(mat)
  res <- matrix(0, n, n, dimnames = list(cn, cn))
  for (i in seq_len(n)) {
    xi <- mat[, i]
    for (j in i:n) {
      d <- JSD(xi, mat[, j])
      res[i, j] <- res[j, i] <- sqrt(d)   # metric化
    }
  }
  as.dist(res)
}

## 거리행렬(필터 후)
dist <- dist.JSD(data, PSEUDO)   # 이하 분석 기본값으로 사용 권장

#####################################################################################

### PAM clustering + CH index over K_RANGE + CH plot

pam.clustering <- function(d, k){
  cluster::pam(d, k, diss=TRUE)$clustering
}

ch_vals <- vapply(K_RANGE, function(k){
  cl <- pam.clustering(dist, k)
  clusterSim::index.G1(t(data), cl, d = dist, centrotypes = "medoids")
}, numeric(1))

k_best <- K_RANGE[which.max(ch_vals)]

pam_best <- cluster::pam(dist, k = k_best, diss = TRUE)
cl_best  <- pam_best$clustering  ##
sil_best <- if (!is.null(pam_best$silinfo)) {
  pam_best$silinfo$avg.width
} else {
  mean(cluster::silhouette(pam_best)[, 3])
}

test <- data.frame(num = K_RANGE, CH = ch_vals)

## 주석 위치(축 범위 비례) 계산
xr <- range(test$num, na.rm = TRUE)
yr <- range(test$CH,  na.rm = TRUE)
x_annot <- xr[2] - 0.25 * diff(xr)   # 더 왼쪽으로: 0.35~0.50로 키우세요
y_annot <- yr[2] - 0.05 * diff(yr)   # 더 아래로: 0.10 등으로 키우세요
sil_lab <- sprintf("Obs. silhouette (k=%d) = %.3f", k_best, sil_best)

ggplot(test, aes(num, CH)) +
  geom_point() +
  geom_segment(aes(xend = num, y = 0, yend = CH,
                   linewidth = ifelse(num == k_best, 1.3, 1.1),
                   color     = ifelse(num == k_best, "orange", "grey")),
               show.legend = FALSE) +
  geom_point(aes(color = ifelse(num == k_best, "orange", "grey"),
                 size  = ifelse(num == k_best, 4, 2.5)),
             show.legend = FALSE) +
  scale_color_identity() + scale_size_identity() + scale_linewidth_identity() +
  labs(x = "Number of clusters", y = "Calinski–Harabasz (CH) index") +
  annotate("text",
           x = max(test$num) - 0.8, y = max(test$CH) * 0.95,
           label = sil_lab,
           size = 4.5, colour = "black") +
  theme(
    plot.background = element_rect(fill='white'),
    rect = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour="black", linewidth=0.5),
    panel.border = element_rect(colour="black", fill=NA, linewidth=0.5),
    plot.title = element_text(size=15, hjust=0.5),
    axis.title = element_text(size=14, color="black", face="bold"), 
    axis.text  = element_text(size=12)
  ) +
  scale_x_continuous(breaks = K_RANGE)

#####################################################################################

### Cluster ↔ Sample mapping (k_best)

cluster_assignment <- data.frame(Subject = colnames(data),
                                 Cluster = cl_best)
cluster_assignment <- cluster_assignment[order(cluster_assignment$Cluster), ]
print(cluster_assignment)

#####################################################################################

### BCA (Between-Class Analysis) with k_best labels

## PCA 입력은 행=샘플 → 전치
obs.pca <- ade4::dudi.pca(t(data), scannf = FALSE, nf = 20)
ng <- length(unique(cl_best))
obs.bet <- ade4::bca(obs.pca, fac = factor(cl_best),
                     scannf = FALSE, nf = max(1, ng - 1))

## 분리 유의성(퍼뮤테이션)
rt <- ade4::randtest(obs.bet, nrepet = 999)
print(rt)  # p-value 확인 (p-value가 작으면(예: <0.05) "그룹 간 구조가 우연 이상으로 존재"함)

## BCA 좌표 시각화
bca_df <- as.data.frame(obs.bet$ls)         # BCA 좌표 (열 개수 = 유효 축 수)
bca_df$Cluster <- factor(cl_best)
n_axes <- ncol(obs.bet$ls)

if (n_axes >= 2) {
  # 2개 이상 축이 있으면 2D 산점도 (축 이름을 표준화해서 사용)
  colnames(bca_df)[1:2] <- c("Axis1","Axis2")
  ggplot(bca_df, aes(Axis1, Axis2, color = Cluster)) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "BCA (between-class analysis)")
} else {
  # 축이 1개뿐이면 1D 분리 플롯
  colnames(bca_df)[1] <- "Axis1"
  ggplot(bca_df, aes(x = Axis1, y = Cluster, color = Cluster)) +
    geom_jitter(height = 0.1, size = 2) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
    theme_minimal() +
    labs(title = "BCA (1D): ng=2 → 축 1개만 존재")
}
## BCA(ng=2 → 축 1개): 그룹 간 분산/그룹 내 분산 비율을 최대화하는 **판별축 1개(Axis1))**만 생성됨. LDA의 1차 축과 유사
# 해석:
#  - Axis1 값이 큰 쪽 = Cluster 2 쪽 성향, 작은 쪽 = Cluster 1 성향.
#  - 점군이 겹치지 않거나 적게 겹치면 군집 분리가 뚜렷.
#  - 0(세로점선)은 중심 기준선일 뿐 절대 임대값이 아님.

#####################################################################################

### PCoA (on filtered JSD) + % variance

pcoa <- ade4::dudi.pco(dist, scannf = FALSE, nf = 3)
scores <- as.data.frame(pcoa$li)
scores$Sample  <- rownames(scores)
scores$Cluster <- factor(cl_best)
pct <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

## 기본 산점도(샘플 라벨)
ggplot(scores, aes(A1, A2, label = Sample, color = Cluster)) +
  geom_point(size = 2) +
  geom_text(vjust = -0.7, size = 2.5, show.legend = FALSE) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_minimal()

###############

# ggscatter 버전(원한다면)
ggscatter(scores, x="A1", y="A2",
          size=1, color="Cluster",
          ellipse=TRUE, ellipse.level=0.8,
          mean.point=TRUE, star.plot=TRUE,
) +
#          label="Sample", font.label=c(9,"plain","black")) +
  xlab(paste0("PC1 (", pct[1], "%)")) +
  ylab(paste0("PC2 (", pct[2], "%)")) +
  theme_minimal()

###############

# (옵션) PCA 고유값 시각화
factoextra::fviz_eig(obs.pca, addlabels = TRUE)

#####################################################################################

### Driver taxa by taxonomy level (Cluster 1 vs 2)

## 0) cl_best ↔ 샘플 정렬 일치시키기
if (is.null(names(cl_best))) {
  names(cl_best) <- colnames(data)              # 이름 없으면 순서대로 부여
} else {
  cl_best <- cl_best[colnames(data)]            # 이름 있으면 컬럼 순서에 맞춰 재정렬
}

## 1) 계통 라벨 파서 (행이름 예: d__Bacteria;p__Firmicutes;...;g__Bacteroides;...)
parse_lineage <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  get <- function(prefix) {
    hit <- grep(paste0("^", prefix, "__"), parts)
    if (length(hit)) sub(paste0("^", prefix, "__"), "", parts[hit[1]]) else NA_character_
  }
  c(kingdom = get("d"), phylum = get("p"), class = get("c"),
    order   = get("o"), family = get("f"), genus = get("g"), species = get("s"))
}
tax_df <- as.data.frame(t(vapply(rownames(data), parse_lineage,
                                 FUN.VALUE = rep(NA_character_, 7))),
                        stringsAsFactors = FALSE, row.names = rownames(data))

## 2) 한 레벨에서 드라이버 상위 n개 계산
top_drivers <- function(level, topn = 10, pseudo = 1e-8) {
  key <- tax_df[[level]]
  key[is.na(key) | key == ""] <- paste0("Unassigned_", level)
  # 같은 레벨 라벨끼리 합산(행 집계)
  agg <- rowsum(data, group = key, reorder = FALSE)
  # 클러스터별 평균
  idx1 <- which(cl_best == 1)
  idx2 <- which(cl_best == 2)
  if (!length(idx1) || !length(idx2)) {
    return(data.frame())                         # 한쪽이 비면 빈 결과
  }
  m1 <- rowMeans(agg[, idx1, drop = FALSE])
  m2 <- rowMeans(agg[, idx2, drop = FALSE])
  diff   <- m2 - m1
  log2fc <- log2((m2 + pseudo) / (m1 + pseudo))
  out <- data.frame(level = rownames(agg),
                    mean_cl1 = m1, mean_cl2 = m2,
                    diff = diff, log2fc = log2fc,
                    row.names = NULL, check.names = FALSE)
  out <- out[order(abs(out$diff), decreasing = TRUE), ]
  head(out, topn)
}

## 3) 레벨별 표 생성 + 콘솔 출력 + CSV 저장
levels_vec <- c("genus","family","order","class","phylum","kingdom","species")
drivers_list <- setNames(vector("list", length(levels_vec)), levels_vec)

out_dir <- getwd()  # 현재 폴더에 저장
for (lev in levels_vec) {
  tab <- top_drivers(lev, topn = 10)
  drivers_list[[lev]] <- tab
  cat("\n=== Top drivers by", toupper(lev), "(|Δ mean| 상위 10) ===\n")
  if (nrow(tab)) print(tab, row.names = FALSE, digits = 4) else cat("No features.\n")
#  fn <- file.path(out_dir, sprintf("drivers_by_%s_top10.csv", lev))
#  utils::write.csv(tab, fn, row.names = FALSE)
#  cat(">> saved:", fn, "\n")
}

# ## (옵션) 하나의 Excel로 저장하고 싶으면 주석 해제 (writexl 필요)
# if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")
# writexl::write_xlsx(
#   c(list(clusters = data.frame(Sample = colnames(data), Cluster = cl_best)),
#     lapply(drivers_list, function(x) { if (is.null(x)) data.frame() else x })),
#   path = file.path(out_dir, "cluster_drivers_by_taxonomy.xlsx")
# )

## 각 csv(drivers_by_<level>_top10.csv)는 해당 레벨에서 클러스터2 - 클러스터1 평균 차이가 큰 순으로 정렬됨.
#  - diff > 0 ⇒ 클러스터 2에서 높은 분류군
#  - diff < 0 ⇒ 클러스터 1에서 높은 분류군
#  log2fc는 방향/크기 참고(0 분모 방지용 pseudo 포함)
