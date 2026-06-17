#install.packages("ggtext")
library(phyloseq)
library(tidyverse)
library(devtools)
library(vegan)
library(glue)
library(ggtext)
library(ggpubr)
library(patchwork)

#### 1. 데이터 변수 입력  ####
setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/no0")

phyloseq = readRDS("physeq_GutLungAxis.rds")  # phyloseq
index = "bray"                  # distance
seed = 42                       # seed값 고정 
type = "Group"              # 비교 group
plot = "PCoA"                   # 차원축소
SampleID = "SampleID"
shap = NULL                     # point 모양 group
type_col = c("#E31A1C", "#1F78B4") # 색 지정 


#### 2. PCoA를 위한 데이터 변환 ####
set.seed(seed)
x.dist <- phyloseq::distance(phyloseq, method = index) # distance계산
ord <- ordinate(phyloseq, plot, index)

# 3. stat - PERMANOVA
myfunc <- function(v1) {
  deparse(substitute(v1))
}
dist <- myfunc(x.dist)
a <- adonis2(as.formula (glue("{dist} ~ {type}")), data=data.frame(sample_data(phyloseq)),
             permutations=9999, method=index)
Perm.p <- a$`Pr(>F)`[1]


#### 4. data 변환 for PCoA  ####
mat <- ord$vectors[, 1:2] %>% as.data.frame() %>%  # PCoA 1, 2 추출
  mutate(SampleID = rownames(.)) %>%
  arrange(SampleID)

meta <- sample_data(phyloseq) %>% data.frame() %>%
  arrange(SampleID) %>%
  mutate(SampleID = rownames(.))

pcoa_df <- inner_join(meta, mat, by = "SampleID") # meta data + PCoA 1,2
## * pcoa_df는 아래에서 계속 사용할 데이터

PC1 <- round(ord$values["Relative_eig"][1,]*100, 1)
PC2 <- round(ord$values["Relative_eig"][2,]*100, 1)


#### 5. PCoA  ####
## main plot
main.plot <- pcoa_df %>% 
  ggplot2::ggplot(aes(x = Axis.1, y=Axis.2)) + 
  ggplot2::geom_vline(xintercept = 0, colour = "grey80") + 
  ggplot2::geom_hline(yintercept = 0, colour = "grey80") + 
  ggplot2::geom_point(aes_string(shape = shap, color=type),alpha = 0.7, size=2.5) +
  lims(x = c(-0.7, 0.9), y = c(-0.6, 0.6)) + 
  ggplot2::stat_ellipse(aes_string(color= type) ) + 
  ggplot2::theme_test() +
  ggplot2::scale_color_manual(values = type_col) +
  ggplot2::labs(
    y =paste0("PCoA2 (", PC2, "%)"),
    x =paste0("PCoA1 (", PC1, "%)")) +
  ggplot2::theme(plot.caption = element_text(hjust = 0)) +
  ggplot2::theme(plot.caption = element_markdown(), 
                 # aspect.ratio=1,
                 legend.position = "bottom", 
                 legend.title = element_blank())+
  theme(legend.position="none",plot.margin=unit(c(0,0,0,0),"points"))+
  ggtext::geom_richtext(label.color = NA, size = 4, fill = NA, 
                        hjust = 0, vjust =0, x = -Inf, y = -Inf, 
                        label= paste0( "**PERMANOVA** *p*-value=", Perm.p ) )
main.plot

ggsave("rarefied_bray1.png", width = 7, height = 6, dpi = 500)

#####################################################################################
#BiocManager::install(c("ggtree", "ggtreeExtra", "treeio", "tidytree"))
library(microbiome)
#devtools::install_github("cafferychen777/MicrobiomeStat")
library(MicrobiomeStat)
library(aplot)

#### 1. import data ####
ps <- readRDS("physeq_GutLungAxis.rds")
data.obj <- mStat_convert_phyloseq_to_data_obj(ps)


#### 2. stat ####
beta.df <- generate_beta_test_single(
  data.obj = data.obj,
  dist.obj = NULL,
  t.level = "2",
  group.var = "Group", 
  dist.name = c('BC', 'Jaccard') 
)
beta.df$p.tab
beta.df$aov.tab

## * 경고메시지(들):
#mStat_calculate_beta_diversity(data.obj, dist.name)에서:
#  It appears the data may not have been rarefied. Please verify.


#### 3. PCoA  ####
data.obj$meta.dat$SubjectID <- rownames(data.obj$meta.dat)
p <- generate_beta_ordination_single(
  data.obj = data.obj,
  dist.obj = NULL,
  pc.obj = NULL,
  subject.var = "SubjectID",
  group.var = "Group",
  strata.var = "SamplingWeek",
  dist.name = c("BC"),
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = FALSE,  #Edit
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5,
)

p

ggsave("rarefied_bray2.png", width = 8, height = 6, dpi = 500)

## * 경고메시지(들):
#mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)에서:
#  It appears the data may not have been rarefied. Please verify.

#####################################################################################

#install.packages("pracma")
library(pracma)

comp_ <- pcoa_df$Group %>% unique() %>% as.vector()
comp_pair <- combn(comp_, 2)

combined_list <- lapply(1:ncol(comp_pair), function(i) comp_pair[, i])

ybox <- 
  ggplot(pcoa_df, aes(x = Group, y = Axis.2)) + 
  geom_boxplot(aes(fill = Group, colour = Group), alpha = 0.3)+
  scale_fill_manual(values = type_col) + 
  scale_color_manual(values = type_col) + 
  lims(y =c(-0.9, 0.9)) +
  stat_compare_means(method = "wilcox.test", tip.length=0.02,
                     label = "p.signif",
                     comparisons = combined_list) + 
  theme_void() + 
  theme(legend.title = element_blank(),
        legend.position = 'none',
        strip.background = element_blank())



xbox <- 
  ggplot(pcoa_df, aes(x = Group, y = Axis.1)) + 
  geom_boxplot(aes(fill = Group, colour = Group), alpha = 0.3)+
  scale_fill_manual(values = type_col) + 
  scale_color_manual(values = type_col) + 
  lims(y =c(-0.9, 0.9)) +
  stat_compare_means(method = "wilcox.test", tip.length=0.02,
                     label = "p.signif",
                     comparisons = combined_list) + 
  theme_void() + 
  theme(legend.title = element_blank(),
        legend.position = 'none',
        strip.background = element_blank())+
  coord_flip()


pcoa <- main.plot %>%
  cowplot::insert_yaxis_grob(ybox, grid::unit(1, "in"), position = "right") %>%
  cowplot::insert_xaxis_grob(xbox, grid::unit(1, "in"), position = "top") %>%
  
  ggdraw()
pcoa

ggsave("rarefied_bray3.png", width = 7, height = 6, dpi = 500)

#####################################################################################

