#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c(
#  "ALDEx2", "DESeq2", "edgeR", "limma", "Maaslin2", "metagenomeSeq", 
#  "SummarizedExperiment", "phyloseq", "biomformat", "lefser"
#))


#install.packages("ggpicrust2")  ## R 4.5 버전에서 이걸로 실행하면 마지막 errorbar 구현 단계에서 오류남. 해결하려면 github로 다운받은 버전 사용.
#remotes::install_github("cafferychen777/ggpicrust2")
#install.packages("ggprism")
#install.packages("GGally")
#install.packages("patchwork")
#install.packages("mutoss")

rm(list = ls())

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(GGally)
library(mutoss)

#####################################################################################

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week2/picrust-MPGA/KO_metagenome_out/pred_metagenome_unstrat.tsv")

metadata <-
  read_delim(
    "../../../../Figures_GutLungAxis/yes0/RUN45_Metadata_rev_GutLungAxis_phyloseq.txt",
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
metadata


colnames(metadata)[1] <- "SampleID"


metadata <- metadata[metadata$Group %in% c("Control", "VNAM"), ]

#####################################################################################

kegg_abundance <-
  ko2kegg_abundance(
    "pred_metagenome_unstrat.tsv"
  )  ## Github 버전 ggpicrust2에선 파일 첫열 이름이 "function."이어야 읽을 수 있음 (점 포함해야됨)
## MPGA 버전 (picrust 2.6 이상)에선 ko:K00001 형태로 나와서 ko: 떼야됨
kegg_abundance

metadata <- metadata[metadata$SampleID %in% colnames(kegg_abundance), ]  ## kegg_abundance에 실제로 존재하는 샘플만 다시 필터
kegg_abundance2 <- kegg_abundance[, metadata$SampleID]

#####################################################################################

metadata$Group <- factor(metadata$Group, levels = c("Control", "VNAM"))

ko_AlDEx <-
  pathway_daa(
    abundance = kegg_abundance2,
    metadata = metadata,
    group = "Group", 
    daa_method = "ALDEx2",
    p.adjust = "none",  ## BKY FDR 하려고 none으로 둠
    select = NULL,
    reference = NULL
  )


ko_AlDEx %>% head()
ko_AlDEx %>% tail()


##**##
## 2) ALDEx2 Wilcoxon 결과만 필터
ko_AlDEx_df <- ko_AlDEx %>%
  dplyr::filter(method == "ALDEx2_Wilcoxon rank test")  ## "ALDEx2_Welch's t test" or "ALDEx2_Wilcoxon rank test"

## 3) BKY(two-stage)로 p-values 보정
p_raw <- ko_AlDEx_df$p_values
bky_res <- mutoss::two.stage(p_raw, alpha = 0.05)  # Benjamini–Krieger–Yekutieli

ko_AlDEx_df$q_bky <- bky_res$adjPValues  ## p_adjust를 찾을 수 없다고 뜨면 이 줄 바꾸면 됨
#ko_AlDEx_df$rej_bky <- bky_res$rejected
ko_AlDEx_df$p_adjust  <- ko_AlDEx_df$q_bky
ko_AlDEx_df$adj_method <- "BKY_two_stage"

## 4) 그룹별 평균으로 log2FC 계산 (VNAM / Control)
ctrl_samples <- metadata$SampleID[metadata$Group == "Control"]
vnam_samples <- metadata$SampleID[metadata$Group == "VNAM"]

ctrl_mean <- rowMeans(kegg_abundance2[, ctrl_samples, drop = FALSE], na.rm = TRUE)
vnam_mean <- rowMeans(kegg_abundance2[, vnam_samples, drop = FALSE], na.rm = TRUE)

## 5) pseudocount 추가
#pseudo <- 1e-6
lfc_vec <- log2( (vnam_mean) / (ctrl_mean) )

lfc_df <- data.frame(
  feature = rownames(kegg_abundance2),
  lfc     = lfc_vec
)


##**##
## 1) sample별 상대풍부도로 변환 (행=pathway, 열=sample 유지)
#rel_mat <- apply(kegg_abundance2, 2, function(x) x / sum(x))

## 2) 그룹별 평균 계산
#ctrl_samples <- metadata$SampleID[metadata$Group == "Control"]
#vnam_samples <- metadata$SampleID[metadata$Group == "VNAM"]

#ctrl_mean_rel <- rowMeans(rel_mat[, ctrl_samples, drop = FALSE], na.rm = TRUE)
#vnam_mean_rel <- rowMeans(rel_mat[, vnam_samples, drop = FALSE], na.rm = TRUE)

## 3) log2FC (VNAM / Control)
#pseudo <- 1e-10
#lfc_rel <- log2((vnam_mean_rel + pseudo) / (ctrl_mean_rel + pseudo))

#lfc_df_rel <- data.frame(
#  feature = rownames(rel_mat),
#  lfc     = lfc_rel
#)
##**##


# 전체 kegg_abundance2로 TSS(relative) 만들기 (subset으로 나누면 안 됨)
rel_all <- sweep(kegg_abundance2, 2, colSums(kegg_abundance2), "/")

ctrl_mean_rel <- rowMeans(rel_all[, ctrl_samples, drop=FALSE], na.rm=TRUE)
vnam_mean_rel <- rowMeans(rel_all[, vnam_samples, drop=FALSE], na.rm=TRUE)

pseudo <- 1e-10
lfc_rel <- log2((vnam_mean_rel + pseudo) / (ctrl_mean_rel + pseudo))

lfc_df_rel <- tibble(feature = rownames(kegg_abundance2), lfc = lfc_rel)



## 4) lfc 열 붙이고, q<0.05 & lfc<-1 로 필터링
ko_AlDEx_df <- ko_AlDEx_df %>%
  dplyr::left_join(lfc_df_rel, by = "feature") %>%
  dplyr::filter(p_adjust < 0.05,
                lfc < 0
                )

ko_AlDEx_df %>% head()
dim(ko_AlDEx_df)
##**##

#####################################################################################

ko_annotation0 <- pathway_annotation(pathway = "KO",
                                     daa_results_df = ko_AlDEx_df,
                                     ko_to_kegg = TRUE)

#####

# KEGG pathway level 1,2,3를 다운 받기
library(KEGGREST)

kegg_abundance2_rev <- kegg_abundance2 %>% rownames_to_column("pathway")
ids <- ko_AlDEx_df$feature

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


#pathway_tab2 <- pathway_tab %>%
#  rownames_to_column("feature") %>%
#  mutate(pathway_class_new = paste(Level1, Level2, sep = "; "))  ## pathway_class로 Level1,2 둘다 가져옴
pathway_tab2 <- pathway_tab %>%
  rownames_to_column("feature") %>%
  mutate(pathway_class_new = Level2)  ## pathway_class로 Level2만 가져옴


# feature 기준으로 ko_annotation에 병합
ko_annotation1 <- ko_annotation0 %>%
  dplyr::select(-pathway_class) %>%
  left_join(
    pathway_tab2[, c("feature", "pathway_class_new", "Level1", "Level2", "Level3")], by = "feature"
#    pathway_tab2[, c("feature", "pathway_class_new")], by = "feature"
    ) %>%
  dplyr::rename(pathway_class = pathway_class_new)


ko_annotation <- ko_annotation1

ko_annotation %>% colnames()


ko_annot_Metabolism <- ko_annotation %>%
  filter(Level1 == "Metabolism") %>%
  dplyr::select(-q_bky, -pathway_name, -pathway_description, -pathway_map, -pathway_class)

write.csv(ko_annot_Metabolism, "../../../picrust2-MPGA_ALDEx2_BKY.csv", row.names = FALSE)

#####################################################################################

ko_annotation <- ko_annotation %>%
  filter(!is.na(pathway_name), !is.na(pathway_class))

Top30 <- ko_annotation %>%
  filter(Level1 == "Metabolism") %>%
  dplyr::filter(lfc < 0) %>%
  arrange(p_adjust) %>%  #Edit: lfc or p_adjust로 정렬. 오름차순으로 정렬됨
#  slice_head(n = 30) %>%
  mutate(pathway_name2 = paste0(pathway_class, " | ", pathway_name))  ##**##

kegg_abundance_t30 <- kegg_abundance2[rownames(kegg_abundance2) %in% Top30$feature, ]

#####################################################################################

p <- pathway_errorbar(abundance = kegg_abundance2,  ## 미리 계산해놓은 lfc 값으로 나타내려고 kegg_abundance_t30 대신 kegg_abundance2를 사용함
                      daa_results_df = Top30,
                      Group = metadata$Group,
                      ko_to_kegg = T,
                      p_values_threshold = 0.05,
                      order = "pathway_class",
                      select = NULL,
                      p_value_bar = T,
                      colors = NULL,
                      x_lab = "pathway_name2",  ##**##
                      pathway_class_position = "none",
                      pvalue_format = "scientific",
                      pvalue_size = 2.4,
                      pvalue_thresholds = c(0.05, 0.01, 0.001, 0.0001),
                      pvalue_star_symbols = c("*", "**", "***", "****"),
                      pathway_class_text_color = "transparent",
                      pathway_class_text_size = 3)

p2 <- p + plot_layout(widths = c(0, 1.4, 0.5, 0.15))  ## 여러 plot을 나타내는거라 그 plot들 순서를 나타냄. 순서대로 전체 plot, 왼쪽부터 차례대로 2-4번째 값을 조절해서 여백 조절.
#p2 <- p2 & theme(axis.text.y = element_text(hjust = 0))
p2 <- p2 & theme(axis.text.y = element_text(hjust = 0, margin = margin(r=5, l=-85)))  ## pathway_class2 글자 왼쪽으로 정렬. margin은 좌우 공백 늘리거나 줄이는 용도.
#p2 <- p2 + plot_annotation(title = "PICRUSt2 - KEGG pathways (Metabolism, Top30)",
#                           theme = theme(plot.title = element_text(hjust = 0.5,
#                                                                   size  = 12,
#                                                                   face  = "bold")))

## Abundance 글자 아래로 내리기
plots <- p2$patches$plots
get_xlab <- function(g){
  if (!inherits(g, "ggplot")) return(NA_character_)
  x <- g$labels$x
  if (is.null(x)) NA_character_ else as.character(x)
}
xlabs <- vapply(plots, get_xlab, character(1))
idx_ab <- which(grepl("Abundance", xlabs, fixed = TRUE))[1]
p2$patches$plots[[idx_ab]] <- p2$patches$plots[[idx_ab]] +
  theme(
    axis.title.x = element_text(margin = margin(t = 20))  ## 숫자 키우면 더 아래로
  )
p2$patches$plots[[3]] <- p2$patches$plots[[3]] +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p2

ggsave("../../../picrust2-MPGA_ALDEx2_BKY_errorbar.png", p2, width = 15, height = 9, dpi = 500)

## 1500/900 png로 저장


### * ko_to_kegg를 FALSE로 두면 오류남
## ** -> 앞서서 pathway_annotation 사용 후, ko_to_kegg를 TRUE로 두는게 오류 날 가능성 적음
### * pathway_annotation을 사용하면 pathway_class 열이 비게 되는데, 이걸 채워줘야 오류 안남
### * Error in ggplot2::annotation_custom(): ! Problem while converting geom to grob. ℹ Error occurred in the 4th layer. Caused by error in UseMethod(): ! 클래스 "c('simpleUnit', 'unit', 'unit_v2')"의 객체에 적용된 'rescale'에 사용할수 있는 메소드가 없습니다 Run rlang::last_trace() to see where the error occurred.
## ** -> update.packages(ask = FALSE) 로 전체 패키지 업데이트 후, remotes::install_github("cafferychen777/ggpicrust2") 로 ggpicrust2 재설치 해서 해결 (25.10.29)

## pvalue_format (default: "smart")
# "numeric": 지수 형태 숫자
# "scientific": 지수 형태 숫자 + 별
# "smart": p-value 범위 (ex. p<0.001) + 별
# "stars_only": 숫자 없이 별만 나타냄
# "combined": 소수점 있는 숫자 형태 + 별

#####################################################################################

