#install.packages("devtools")
#devtools::install_github("grunwaldlab/metacoder")
library(phyloseq)
library(metacoder)
library(ggplot2)
library(dplyr)

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/no0")

ps <- readRDS("physeq_GutLungAxis.rds") %>%
  transform_sample_counts(function(x) x/sum(x) ) %>% 
  subset_samples(Group %in% c("Control", "VNAM"))
ps # 770 taxa and 17 samples
## * 그룹 3개 이상 지정 가능 (오류 안나긴 함..)

# read가 0인 샘플, taxa제거
ps.f = prune_samples(sample_sums(ps)>0, ps) 
#ps.f = prune_taxa(rowSums(otu_table(t(ps.f))) > 0, ps.f)  # 에러 -> 지금까지 모든 데이터에서 오류 났었음 (25.11.02)

#####################################################################################

obj <- parse_phyloseq(ps.f)

obj %>%
  filter_taxa(taxon_ranks == "Order", supertaxa = TRUE)

#####################################################################################

obj$data$tax_data <- calc_obs_props(obj, "tax_data") # 각 taxa 계산
obj$data$tax_abund <- calc_taxon_abund(obj, "otu_table") # 샘플의 abundance 계산
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", 
                                   groups = obj$data$sample_data$Group,  ##
                                   cols = obj$data$sample_data$sample_id)

set.seed(42)
obj %>% 
  filter_taxa(taxon_ranks == "Order", supertaxa = TRUE)

#####################################################################################

obj %>%
  filter_taxa(taxon_ranks == "Order", supertaxa = TRUE)

#####################################################################################

obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id,
                                      groups = obj$data$sample_data$Group)  ##
print(obj$data$diff_table)

obj %>% 
  filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% 
  heat_tree(node_label = taxon_names,
            node_size_axis_label = "Number of OTUs",
            node_size = n_obs,
            node_color_axis_label = "Mean difference",
            node_color = mean_diff,
            node_color_range= c("#E31A1C","grey90",  "#1F78B4"), 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford",
            background_color = "white",                 # ← 흰배경
            output_file = "tree1_genus_num of otus.png")

## * 경고메시지(들):
#1: There is no "taxon_id" column in the data set "3", so there are no taxon IDs. 
#2: The data set "4" is named, but not named by taxon ids.

#####################################################################################

ps <- readRDS("physeq_GutLungAxis.rds") %>%  ##
  subset_samples(Group %in% c("Control", "VNAM"))  ##
ps.f = prune_samples(sample_sums(ps)>0, ps) 
obj2 = parse_phyloseq(ps.f)

## * The following 641 of 765 (83.8%) input indexes have `NA` in their classifications:
#1, 2, 3, 4, 5, 6, 7, 8, 9 ... 757, 758, 759, 760, 761, 762, 763, 764, 765


obj2$data$tax_data <- calc_obs_props(obj2, "tax_data") # 각 taxa 계산
obj2$data$tax_abund <- calc_taxon_abund(obj2, "otu_table") # 샘플의 abundance 계산
obj2$data$tax_occ <- calc_n_samples(obj2, "tax_abund", 
                                    groups = obj2$data$sample_data$Group,  ##
                                    cols = obj$data$sample_data$sample_id)
obj2$data$tax_data
obj2$data$tax_abund
obj2$data$tax_occ

## * No `cols` specified and no numeric columns can be found.
## * 경고메시지(들):
#  No cols specified. No calculation will be done.

obj2$data$diff_table <- calc_diff_abund_deseq2(obj2, data = "tax_abund",
                                               cols = obj$data$sample_data$sample_id,
                                               groups =  obj$data$sample_data$Group)  ##

print(obj2$data$diff_table)


## * 경고메시지(들):
#DESeqDataSet(se, design = design, ignoreRank)에서:
#  some variables in design formula are characters, converting to factors

obj2 %>% 
  filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color_axis_label = "Log2 fold change",
            node_size_axis_label = "Number of OTUs",
            node_color = ifelse(is.na(padj) | padj > 0.05, 0, log2FoldChange),
            node_color_range= c("#E31A1C","grey90",  "#1F78B4"), 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford",
            background_color = "white",
            output_file = "tree2_genus_lfc.png")

## * 경고메시지(들):
#1: There is no "taxon_id" column in the data set "3", so there are no taxon IDs. 
#2: The data set "4" is named, but not named by taxon ids

#####################################################################################

