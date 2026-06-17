#BiocManager::install("microbiome")
library(phyloseq)
library(ggplot2)
library(microbiome)
library(viridis)
library(RColorBrewer)

setwd("C:/Users/tr432/Downloads/filter_GutLungAxis_Week2")

ps <- readRDS("physeq_GutLungAxis.rds")
ps_genus <- tax_glom(ps, taxrank = "Genus")

ps_genus.rel <- microbiome::transform(ps_genus, "compositional")
ps_genus.rel2 <- prune_taxa(taxa_sums(ps_genus.rel) > 0, ps_genus.rel)


## y축 라벨을 Genus 이름으로 바꾸기 (기본은 OTU ID)
tax <- as.data.frame(tax_table(ps_genus.rel2))
genus_names <- as.character(tax$Genus)
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified_Genus"
taxa_names(ps_genus.rel2) <- genus_names


core.taxa.standard <- core_members(ps_genus.rel2, detection = 1/1000, prevalence = 50/100)
head(core.taxa.standard)


pseq.core <- core(ps_genus.rel, detection = 0.1/100, prevalence = 50/100)
pseq.core


prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)


p.core <- plot_core(ps_genus.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) 
p.core

ggsave("core_microbiome_GutLungAxis.png", width = 9, height = 7, dpi = 500)


taxonomy <- as.data.frame(tax_table(ps_genus.rel2))


core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)
head(core_taxa_id$Family)

write.csv(core_taxa_id, file = "core_taxa_id_GutLungAxis.csv", row.names = TRUE)
