#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}

#BiocManager::install("microbiomeMarker")

rm(list = ls())

library(microbiomeMarker)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(phyloseq)

######################################################################################

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week2")

ps <- readRDS("physeq_GutLungAxis.rds")
ps

### (선택) 그룹 필터링
#ps_fru <- subset_samples(ps, Group %in% c("24h_control","Fructose_GCU"))
#ps_fru
#table(sample_data(ps_fru)$Group)

######################################################################################

tax <- as.data.frame(tax_table(ps))

for (i in seq_len(ncol(tax))) {
  rank_name <- colnames(tax)[i]
  tax[[i]] <- as.character(tax[[i]])
  tax[[i]][is.na(tax[[i]]) | tax[[i]] == "" | tax[[i]] == "NA"] <- paste0(rank_name, "_unclassified")
}

tax_table(ps) <- as.matrix(tax)

######################################################################################

set.seed(0)

mm_test <- run_lefse(ps,
#                     norm = 'CPM', 
                     group = 'Group',
                     taxa_rank = 'all',  #Edit: plot_ef_bar 그릴땐 'Genus'로, plot_cladogram 그릴땐 'all'로
                     kw_cutoff = 0.001,
                     wilcoxon_cutoff = 0.001,
                     multigrp_strat = TRUE,
                     lda_cutoff = 4)

marker_table(mm_test)

######################################################################################

plot_ef_bar(mm_test)

######################################################################################
## marker 개수 확인
mt <- as.data.frame(marker_table(mm_test))

nrow(mt)
table(mt$enrich_group)
summary(mt$ef_lda)


plot_cladogram(
  mm_test,
  color = c('Control' = 'darkgreen', 'VNAM' = 'red'),
  only_marker = TRUE,
  clade_label_level = 3,
  clade_label_font_size = 2.5
)

## * mm_test$marker_table$enrich_group에 있는 두 그룹명 그대로 써야됨

######################################################################################

#BiocManager::install("microbial")
#BiocManager::install("microbiome")

library(microbial) 
library(microbiome)


# lefse 계산
res <- ldamarker(ps, group="Group")  ##

# plot 그리기
lefse_plot <- plotLDA(res,
                      group=c("Control","VNAM"),
                      lda=5,
                      pvalue=0.05,
                      padj = "BH",
                      color = c("red","green"),
                      fontsize.x = 7,
                      fontsize.y = 8)
lefse_plot

######################################################################################

