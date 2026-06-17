rm(list = ls())

#install.packages("ggstar")
#install.packages("TDbook")
#install.packages("reshape")
library(phyloseq)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)
library(reshape)


setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week0")

ps <- readRDS("physeq_GutLungAxis.rds")
ps


## read count to relative abundance
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x) )


##**##
tt <- as.data.frame(tax_table(ps.rel))

# 모든 rank 컬럼에서 "k__/p__/c__/o__/f__/g__/s__" 접두어 제거
tt2 <- tt %>%
  mutate(across(everything(), ~ str_remove(as.character(.x), "^[a-z]__")))

tax_table(ps.rel) <- tax_table(as.matrix(tt2))
##**##


## Bacteria만
ps.rel <- subset_taxa(ps.rel, Kingdom == "Bacteria")

## Genus까지로 collapse (Genus가 같은 taxa들을 합산)
ps.rel <- tax_glom(ps.rel, taxrank = "Genus", NArm = TRUE)

## (선택) Unassigned/빈 Genus 제거
ttg <- as.data.frame(tax_table(ps.rel))
keep <- !is.na(ttg$Genus) & ttg$Genus != "" & ttg$Genus != "Unassigned"
ps.rel <- prune_taxa(keep, ps.rel)


##**##
tt <- as.data.frame(tax_table(ps.rel))

# 금지 키워드 설정
keep_taxa <- rownames(tt)[
  !str_detect(tt$Family, "^(CAG|UBA)") &
    !str_detect(tt$Genus, "^(UBA|HUN|CAG|COE|C-53|MD|TWA)")
  ]  #Edit: 직접 plot 보고 판단해서 금지 키워드 설정

ps.rel <- prune_taxa(keep_taxa, ps.rel)
ps.rel
##**##

## 데이터 일부만 추출 ##
ps.rel.bac <- subset_taxa(ps.rel, Kingdom == "Bacteria" & 
                            Family != "Unclassified Bacteria")  #Edit
#myTaxa = names(sort(taxa_sums(ps.rel.bac), decreasing = TRUE))
myTaxa = names(sort(taxa_sums(ps.rel.bac), decreasing = TRUE)[1:50])  #Edit
ps.50 = prune_taxa(myTaxa, ps.rel.bac)
ps.50

#####################################################################################
## 계통수 그리기
# 기본적인 계통수 그리기. extrenal node의 색은 Family별로 구분. 이때 기본적인 OTU단위로 node가 생성됨

ggtree(ps.50, layout="circular", open.angle=10, ) +
  geom_tippoint(mapping=aes(color=Family), size=1.5, show.legend=T)

#####################################################################################
# 추가적으로 각 OTU의 relative abundance의 정보를 boxplot으로 표현

p <- ggtree(ps.50, layout="circular", open.angle=10, ) +
  geom_tippoint(mapping=aes(color=Family), size=1.5, show.legend=F) +  #Edit
  geom_fruit(geom=geom_boxplot,
             mapping = aes(
               y=OTU,
               x=Abundance,
               fill=Family),  #Edit
             size=.2,
             outlier.size=0.5,
             outlier.stroke=0.5,
             outlier.shape=21,
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3),
             grid.params=list())

#p

# otu가 각 그룹에서 얼마나 풍부하게 존재하는지 알아보기 위해 타일 형태의 그림(geom_tile)을 추가.이때 ggnewscale::new_scale_fill()을 적은 후 그 밑에 그리고 싶은 내용을 추가

p2 <- p + new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(y=OTU,
                         x = `Group`,  #Edit
                         fill=`Group`,  #Edit
                         alpha = Abundance),
             color = "grey30", offset = 0.04, size = 0.02) +
  scale_alpha_continuous(range=c(0, 1),guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))

#p2

# 각 node에 OTU ID나 Phylum정보를 붙이고 싶다면 geom_tiplab(aes(label= {otu혹은 taxa level})) 을 사용

p2 + geom_tiplab(aes(label=Genus), size = 2.5)  #Edit

ggsave("ggtree1_top50.tiff", width = 15, height = 8, dpi = 800)

#####################################################################################
# 아래처럼 트리의 구조를 바꾸어도 동일한 내용을 출력할 수 있음

p3 <- ggtree(ps.50, layout="rect", open.angle=10, ) +
  geom_tippoint(mapping=aes(color=Family), size=1.5, show.legend=F) +  #Edit
  geom_fruit(geom=geom_boxplot,
             mapping = aes(
               y=OTU,
               x=Abundance,
               fill=Family),  #Edit
             offset = 0.35,  #Edit: plot 간 여백 설정
             pwidth = 1,  #Edit: plot 폭 설정
             size=.2,
             outlier.size=0.5,
             outlier.stroke=0.5,
             outlier.shape=21,
             axis.params=list(
               axis       = "x",
               text.size  = 1.8,
               hjust      = 1,
               vjust      = 0.5,
               nbreak     = 3),
             grid.params=list()) +
  new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(y=OTU,
                         x = `Group`,  #Edit
                         fill=`Group`,  #Edit
                         alpha = Abundance),
             color = "grey30", offset = 0.04, size = 0.02) +
  scale_alpha_continuous(range=c(0, 1),guide=guide_legend(keywidth = 0.3, keyheight = 0.3, order=5))

#p3

p3 + geom_tiplab(aes(label=Genus), size = 2.5)  #Edit

ggsave("ggtree2_top50.tiff", width = 15, height = 8, dpi = 800)

#####################################################################################

