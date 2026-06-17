#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library(qiime2R)

rm(list = ls())

#####################################################################################
### Reading Artifacts (.qza)

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")

SVs <- read_qza("filtered_table_freq-min-1000.qza")

names(SVs)

SVs$data[1:5,1:5]

SVs$type
# [1] "FeatureTable[Frequency]"

#####################################################################################
### Reading Metadata

metadata <- read_q2metadata("RUN45_Metadata_rev_GutLungAxis_phyloseq.txt")
## * metadata 파일 SampleID 행 아래에 #q2:types (categorical/numeric) 추가해야됨
# ** numerical 열에 N/A 같은 문자가 있으면 오류남. 다 지워야 됨 (숫자로 바꾸거나)
head(metadata)

#####################################################################################
### Reading Taxonomy

taxonomy <- read_qza("gg2_taxonomy.qza")
head(taxonomy$data)

taxonomy <- parse_taxonomy(taxonomy$data)
head(taxonomy)

#####################################################################################
### Creating a Phyloseq Object

physeq <- qza_to_phyloseq(
  features="filtered_table_freq-min-1000.qza",
  tree="rooted_tree.qza",
  taxonomy = "gg2_taxonomy.qza",
  metadata = "RUN45_Metadata_rev_GutLungAxis_phyloseq.txt"
)
physeq


saveRDS(physeq, file="physeq_GutLungAxis_freq-min-1000.rds")

#####################################################################################
#####################################################################################
### Visualizations
## Plotting a Heatmap

rm(list = ls())

library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("RUN45_Metadata_rev_GutLungAxis_phyloseq.txt")
SVs <- read_qza("filtered_table_freq-min-1000.qza")$data
taxonomy <- read_qza("gg2_taxonomy.qza")$data

SVs <- apply(SVs, 2, function(x) x/sum(x)*100)  ## 상대적 풍부도로 변환


## SVs + taxonomy 붙인 뒤, Genus 추출
sv_long <- SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  left_join(taxonomy, by = c("Feature.ID" = "Feature.ID")) %>%
  mutate(
    Genus_raw = str_extract(Taxon, "g__[^;]+"),
    Genus     = str_replace(Genus_raw, "g__", ""),
    Genus     = if_else(is.na(Genus) | Genus == "", "Unassigned", Genus)
  ) %>%
  filter(!is.na(Genus), Genus != "NA") %>%
  tidyr::gather(-Feature.ID, -Taxon, -Genus,
                key = "SampleID", value = "Abundance") %>%
  filter(!is.na(Abundance),
         str_detect(Abundance, "^[0-9\\.eE+-]+$")) %>%  ## 숫자인 것만
  mutate(
    Abundance = as.numeric(Abundance)
  )


## Genus별 Toal Abundance 계산 -> 상위 30개 Genus
top_genus <- sv_long %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 30) %>%  ## 상위 30개 Genus
  pull(Genus)
top_genus

sv_genus_plot <- sv_long %>%
  group_by(SampleID, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate( Feature = if_else(Genus %in% top_genus, Genus, "Remainder") ) %>%
  left_join(metadata, by = c("SampleID" = "SampleID")) %>%
  filter(!is.na(Group)) %>%  #Edit
  mutate(NormAbundance = log10(Abundance + 0.01))  ## pseudocount


## Feature(Genus)를 total abundance 기준으로 오름차순 정렬
order_df <- sv_genus_plot %>%
  group_by(Feature) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(TotalAbundance)  ## 내림차순으로 바꾸려면 desc() 씌우면 됨

## Unassigned 제거
order_levels <- setdiff(order_df$Feature, "Unassigned")

## Remainder를 맨 아래(= levels의 첫 번째)로
order_levels <- c("Remainder",
                  setdiff(order_levels, "Remainder"))

sv_genus_plot$Feature <- factor(sv_genus_plot$Feature,
                                levels = order_levels)


sv_genus_plot <- sv_genus_plot %>% 
  filter(!is.na(Feature), Feature != "NA", Feature != "Unassigned")


## Heatmap (y축 = Genus (+ Remainder))
sv_genus_plot %>%
  ggplot(aes(x = SampleID, y = Feature, fill = NormAbundance)) +
  geom_tile() +
  facet_grid(~ Group, scales = "free_x") +  #Edit
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(name = "log10(% Abundance)")

ggsave("heatmap_genus_freq-min-1000.png", width = 9, height = 6, dpi = 500)

#####################################################################################

