rm(list = ls(all.names = TRUE))
gc()
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F)

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

setwd('C:/Users/tr432/Downloads/filter_GutLungAxis_Week2/picrust-MPGA/KO_metagenome_out/pred_metagenome_unstrat.tsv')

## ggpicrust2의 패키지를 이용해서 KO abundance를 kegg pathway abundance로 변환
kegg_abundance <- read.csv("pred_metagenome_unstrat.tsv", sep = '\t', check.names = FALSE)  ## Github 버전 ggpicrust2에선 row.names = 1 빼야됨 <25.11.30>
head(kegg_abundance)
dim(kegg_abundance)

## (선택) ko가 아닌 K 형태일 때 사용
#kegg_abundance <- rownames_to_column(kegg_abundance)


kegg_abundance <- ko2kegg_abundance(data = kegg_abundance)
kegg_abundance

#####################################################################################

# kegg_abundance에서도 원하는 샘플만 추출
metadata0 <- read.csv("../../../../Figures_GutLungAxis/yes0/RUN45_Metadata_GutLungAxis.txt", sep = '\t', check.names = FALSE)
colnames(metadata0)[1] <- "SampleID"
metadata <- metadata0[metadata0$Group %in% c("Control", "VNAM"), ]
metadata <- metadata[metadata$SampleID %in% colnames(kegg_abundance), ]  ## df에 실제로 존재하는 샘플만 다시 필터
kegg_abundance2 <- kegg_abundance[, metadata$SampleID]
dim(kegg_abundance2)

#####################################################################################

Group <- factor(metadata$Group)
Level <- levels(Group)

## ALDEx2 분석 ##
# 1) Aldex2분석을 위해 OTU read count 데이터를 CLR로 normalization
class(metadata$Group)
ALDEx2_object <- ALDEx2::aldex.clr(round(kegg_abundance2), as.character(metadata$Group))
# clr 변환은 각 샘플의 모든 taxa값이 0이며, taxa간의 상대적 비율 값이 그대로 보존된다

# 2) Aldex 통계 분석
ALDEx2_results <- ALDEx2::aldex.ttest(ALDEx2_object, paired.test = FALSE, verbose = FALSE)
# t.test결과와 wilcoxon rank sum test결과를 반환한다 
# 각 taxa가 그룹내에서 유의한 수준으로 차이가 '있는지 없는지' 판별

# 3) Effect size계산
ALDEx2_effect <- ALDEx2::aldex.effect(ALDEx2_object)
# 각 taxa가 차이가 있다면(유의하다면), '얼마나' 차이가 있는지 판별


# 4) 결과물 정리 
p_values_df <- data.frame(
  feature = rep(rownames(ALDEx2_results), 2), 
  method = c(rep("ALDEx2_Welch's t test", nrow(ALDEx2_results)), 
             rep("ALDEx2_Wilcoxon rank test", nrow(ALDEx2_results))),
  group1 = rep(Level[1], 2 * nrow(ALDEx2_results)),
  group2 = rep(Level[2], 2 * nrow(ALDEx2_results)),
  p_values = c(ALDEx2_results$we.ep, ALDEx2_results$wi.ep),
  effect = ALDEx2_effect$effect)

# (선택) 5-1) P- value값 보정 (fdr 사용시)
# BH사용시, 기존 p_values_df에 wi.eBH를 사용 
#adjusted_p_values <- data.frame(
#  feature = p_values_df$feature,
#  p_adjust = p.adjust(p_values_df$p_values, method = "fdr"))


##**##
## (선택) 5-2) P-value 보정: BKY(two-stage Benjamini–Krieger–Yekutieli)
# Wilcoxon만 추출
p_wi <- p_values_df %>%
  filter(method == "ALDEx2_Wilcoxon rank test")

bky_wi <- mutoss::two.stage(pValues = p_wi$p_values,
                             alpha   = 0.05)

p_wi$q_bky <- bky_wi$adjPValues

ko_AlDEx_result <- p_wi %>%
  dplyr::rename(p_adjust = q_bky)
##**##

# (선택) 6) 결과물 최종 정리 (BKY 사용 안할 때)
#ko_AlDEx_result <- cbind(p_values_df, p_adjust = adjusted_p_values$p_adjust)
#ko_AlDEx_result

#####################################################################################

# "ALDEx2_Wilcoxon rank test"결과 중에서 p-value값 보정이 유의한 것만 추출
ko_AlDEx_result2 <- ko_AlDEx_result[ko_AlDEx_result$method == "ALDEx2_Wilcoxon rank test" & ko_AlDEx_result$p_adjust < 0.05 & ko_AlDEx_result$effect < 0, ]
dim(ko_AlDEx_result2)


# KEGG pathway level 1,2,3를 다운 받기
library(KEGGREST)

kegg_abundance2_rev <- kegg_abundance2 %>% rownames_to_column("pathway")
ids <- kegg_abundance2_rev$pathway

# kegg database 다운로드
Kegg_results <- list()
for ( i in ids) {
  Kegg_results[[i]] <- tryCatch(keggGet(i), error=function(e) NULL) # https://www.biostars.org/p/366463/
}

# level 1,2,3를 담을 data.frame만들기
keg_ids <- names(Kegg_results)
pathway_tab <-  data.frame(row.names = keg_ids)

# for 문을 이용해서 level데이터 가져오기
for (i in keg_ids){
  pathway_tab[i, "Level1"] <- strsplit( Kegg_results[[i]][[1]]$CLASS, "; ")[[1]][1]
  pathway_tab[i, "Level2"] <- strsplit( Kegg_results[[i]][[1]]$CLASS, "; ")[[1]][2]
  pathway_tab[i, "Level3"] <- Kegg_results[[i]][[1]]$PATHWAY_MAP
  
}


pathway_tab2 <- pathway_tab %>% rownames_to_column("feature")
ko_ann <- merge(ko_AlDEx_result2, pathway_tab2, by = "feature")

#####################################################################################

# Pathway 1,2,3을 하나의 이름으로 표시하기
sum(is.na(ko_ann$Level3))
pathway_table <- ko_ann %>%
  filter(!is.na(Level3)) %>%
  na.omit %>%  # 44
#  mutate(Level_123 = paste0(Level1, "; ", Level2," - ", Level3))
  mutate(Level_123 = paste0(Level2," - ", Level3))
rownames(pathway_table) <- pathway_table$feature


# annotation이 완료된 kegg만 추출하기
kegg_abundance3 <- kegg_abundance2[pathway_table$feature, ]


# kegg abundance를 relative abundance로 바꾸기
#col_sums <- colSums(kegg_abundance3)
#relative_abundance <- t(t(kegg_abundance3) / col_sums)

# heatmap을 위한 데이터 생성
#table <- merge(pathway_table, relative_abundance, by = "row.names")
#rownames(table) <- table$Level_123


# (추천) 원값 + log 변환
mat_hm <- log1p(as.matrix(kegg_abundance3))  # 또는 as.matrix(kegg_abundance3)

# heatmap을 위한 데이터 생성 (원값 사용)
table <- merge(pathway_table, mat_hm, by = "row.names")
rownames(table) <- table$Level_123

#####################################################################################

## Top 30 선택
table <- table %>%
  filter(!is.na(Level2), !is.na(Level3))

table <- table %>%
  filter(Level1 == "Metabolism") %>%
  arrange(effect)# %>%  #Edit: 원하는 기준으로 정렬 (ex. p_adjust, effect)
#  slice_head(n = 30)  #Edit: 원하는만큼 Top 개수 선택

#####################################################################################

# heatmap 그리기
library(pheatmap)
#p <- pheatmap(table[, 13:33])  #Edit: [, 13:(table 열 개수만큼)] 설정
p <- pheatmap(table[, metadata$SampleID, drop = FALSE])

#####################################################################################

annotation_row <- data.frame(
  row.names  = rownames(table),
  p_adjust = table$p_adjust,
  Effect = table$effect
)

metadata$Group %>% table
annotation_row <- annotation_row[order(annotation_row$Effect, decreasing = F), ,
                                 drop = FALSE # rownames 사라지는걸 막음
]

#####################################################################################

annotation_col <- data.frame(
  row.names = metadata$SampleID,
  Group =  metadata$Group  #Edit
)

annotation_col$Group <- factor(annotation_col$Group , levels = c("Control", "VNAM") )

annotation_col <- annotation_col[order(annotation_col$Group, decreasing = FALSE), , drop = FALSE  ]

#####################################################################################

ann_colors = list(
  p_adjust = colorRampPalette(c("white", "green"))(100),
#  Effect = colorRampPalette(c("#FFC0CB", "white", "#8CD0EC"))(100),
  Effect = colorRampPalette(c("#FFC0CB", "white"))(100),
  Group = c("Control" = 'firebrick',"VNAM" = 'royalblue')
)

#####################################################################################

#df <- table[, 13:33]  #Edit: [, 13:(table 열 개수만큼)] 설정
df <- table[, metadata$SampleID, drop = FALSE]
colnames(df)
df <- df[rownames(annotation_row), rownames(annotation_col)]

sum(annotation_col$Group == "Control")  # gaps_col에 입력
sum(table$effect < 0)                  # gaps_row에 입력


plot <- pheatmap::pheatmap(mat = as.matrix(df),
                           color = colorRampPalette(c("blue", "white", "red"))(100),
                           annotation_col = annotation_col, 
                           annotation_row = annotation_row, 
                           annotation_colors = ann_colors,
                           scale = "row", # 행별로 정규화
                           cluster_rows = F,  
                           cluster_cols = F,
                           gaps_col = 6, # n번째 열에서 heatmap 분리
                           ## 두 그룹 중 첫번째 그룹 샘플개수만큼 입력 (annotation_col 데이터 확인)
                           gaps_row = 28,# n번째 행에서 heatmap 분리
                           ## table 데이터 확인해서 effect가 0보다 작은 샘플개수만큼 입력
                           legend = T,
                           border_color=NA)

png("../../../picrust2-MPGA_ALDEx2_BKY_heatmap.png", width = 18, height = 9, units = "in", res = 500)
plot
dev.off()

# 1600 / 800


# 순서 정렬
pathway_table2 <- merge(pathway_table, table,
                        by.x = "feature",    # pathway_table 에 있는 열 이름
                        by.y = "feature",  # table 에 있는 열 이름
                        all = FALSE          # 겹치는 것만 (inner join)
                        )
pathway_table2 <- pathway_table2[, !grepl("\\.y$", names(pathway_table2))]
names(pathway_table2) <- sub("\\.x$", "", names(pathway_table2))

lev = pathway_table[order(pathway_table2$effect), "Level3"]
pathway_table2$Level_123 <- factor(pathway_table2$Level3, level = lev)
# ggplot
pathway_table2 %>% 
  mutate(Group = if_else(effect<0, "Control", "VNAM")) %>% 
  filter(abs(effect)>0.5) %>%
  ggplot(aes(x = effect, y = Level_123, fill = Group)) + 
  geom_col() + theme_classic() + 
  labs(y = NULL, x = 'Effect size')+
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by=0.5))

# width 1440 / height 830

#####################################################################################

