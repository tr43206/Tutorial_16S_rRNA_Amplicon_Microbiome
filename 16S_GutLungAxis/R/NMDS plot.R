rm(list = ls())

library(vegan)
library(phyloseq)

setwd("C:/Users/tr432/Downloads")

ps <- readRDS("ps.ng.tax.rds")
ps
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4710 taxa and 474 samples ]
# sample_data() Sample Data:       [ 474 samples by 30 sample variables ]
# tax_table()   Taxonomy Table:    [ 4710 taxa by 6 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4710 tips and 4709 internal nodes ]


otu_table(ps)[1, ]

#####################################################################################

# 데이터 크기를 줄여서 진행
ps_rs = prune_taxa(taxa_sums(ps) > 100, ps)
ps_prop <- transform_sample_counts(ps_rs, function(x) 100*x/sum(x))  ## 상대적 풍부도로 바꿈

#####################################################################################
## Sample별 NMDS plot 그리기

set.seed(123)  ## 25.11.03 추가
otu <- otu_table(ps_prop)
ord <- ordinate(physeq = ps_prop,
                method = "NMDS",  # Default is "DCA" or "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
                distance = "bray")  # Default is "bray" or "unifrac", "wunifrac", "jaccard",
p <- plot_ordination(ps_prop, ord, "samples", color="sample_type")
p

#####################################################################################
## ggplot으로 그리기

samples <- data.frame(scores(ord$points))
samples$SampleType<- get_variable(ps_prop, "sample_type")  ## sampledata에 "sample_type"이라는 열을 가져와서 열 추가 

head(samples)
#                                          MDS1        MDS2 SkinType
#1927.SRS014341.SRX020546.SRR043648 -0.81474818  0.16662040  Vaginal
#1927.SRS014837.SRX020546.SRR043659 -1.05460391 -0.33291868  Vaginal
#1927.SRS015389.SRX020546.SRR043661  0.57928105 -0.05124357    Stool
#1927.SRS014627.SRX020546.SRR043662 -0.92117159 -0.28196423  Vaginal
#1927.SRS014923.SRX020546.SRR043675  0.47449253 -0.16889302    Stool
#1927.SRS015127.SRX020546.SRR043678 -0.05396112  0.34111887     Oral

## * 결과값이 다름. -> (?) 시드값 설정할때마다 달라지나..?


library(ggplot2)
p <- ggplot(samples) +
  geom_point(aes(x = MDS1, y = MDS2, color=SampleType)) +
  theme_bw()
p

#####################################################################################

ggplot(samples) +  # sets up the plot. brackets around the entire thing to make it draw automatically
  geom_point(aes(x = MDS1, y = MDS2, color=SampleType)) +  # puts the site points in from the ordination, shape determined by site, size refers to size of point
  stat_ellipse(level = 0.9, aes(x = MDS1, y = MDS2, fill=SampleType),geom = "polygon", alpha = 0.3)+
  theme_bw()+
  theme(legend.position = "none")

#####################################################################################
## phylum 별로 NMDS 분석 plot그리기

spps <- as.data.frame(scores(ord, choices=c(1,2),display=c("species")))
spps_species <- as.data.frame(tax_table(ps_prop))

head(spps_species)
#             Domain        Phylum       Class         Order         Family       Genus
#9410494158 Bacteria Bacteroidetes Bacteroidia Bacteroidales Prevotellaceae  Prevotella
#9410494576 Bacteria Bacteroidetes Bacteroidia Bacteroidales Prevotellaceae  Prevotella
#9410491420 Bacteria Bacteroidetes Bacteroidia Bacteroidales Bacteroidaceae Bacteroides
#9410491569 Bacteria Bacteroidetes Bacteroidia Bacteroidales Bacteroidaceae Bacteroides
#9410491595 Bacteria Bacteroidetes Bacteroidia Bacteroidales Bacteroidaceae Bacteroides
#9410491783 Bacteria Bacteroidetes Bacteroidia Bacteroidales Bacteroidaceae Bacteroides

spps$Phylum <- factor(spps_species$Phylum)  # making a column with species names

head(spps)
#               NMDS1       NMDS2        Phylum
#9410491526 0.6314447 -0.08904021 Bacteroidetes
#9410491516 0.6231495 -0.10382150 Bacteroidetes
#9410492612 0.4698253 -0.37411464 Bacteroidetes
#9410491521 0.6435790 -0.04185019 Bacteroidetes
#9410491824 0.5839583 -0.02245883 Bacteroidetes
#9410491817 0.5585791 -0.16516024 Bacteroidetes

#####################################################################################

p <- ggplot() +
  geom_point(data=spps, aes(y=NMDS2, x=NMDS1, color=Phylum)) +
  theme_bw()
p

#####################################################################################
## phylum + Sample type별 NMDS plot그리기

p <- ggplot(samples) +  # sets up the plot. brackets around the entire thing to make it draw automatically
  geom_point(aes(x = MDS1, y = MDS2, shape=SampleType)) +  # puts the site points in from the ordination, shape determined by site, size refers to size of point
  geom_point(data=spps, aes(y=NMDS2, x=NMDS1, color=Phylum)) +
  theme_bw()
p

#####################################################################################
## SampleType에 ellipse 추가

p <- ggplot(samples) +
  geom_point(aes(x = MDS1, y = MDS2, shape=SampleType)) +
  geom_point(data=spps, aes(y=NMDS2, x=NMDS1, color=Phylum)) +
  stat_ellipse(level = 0.9, aes(x = MDS1, y = MDS2, fill=SampleType),geom = "polygon", alpha = 0.3) +
  theme_bw()
p

#####################################################################################

