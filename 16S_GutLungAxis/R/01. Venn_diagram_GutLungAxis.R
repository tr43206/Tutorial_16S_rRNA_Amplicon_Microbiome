library(phyloseq)
library(dplyr)

rm(list = ls())

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")

ps <- readRDS("physeq_GutLungAxis_freq-min-1000.rds")

#####################################################################################

meta <- data.frame(sample_data(ps))
otu <- data.frame(otu_table(ps))

colnames(otu) <- gsub("\\.", "-", colnames(otu))
common <- intersect(rownames(meta), colnames(otu))
meta <- meta[common, , drop = FALSE]
otu  <- otu[, common, drop = FALSE]

table <- merge(meta, t(otu), by = "row.names")

#####################################################################################

Control <- table %>%
  filter(Group %in%  "Control") %>%
    dplyr::select(which(colSums(.[28:5416]) > 0)) %>%  # feature table이 column 6~643에 위치, 이 중 전체 합이 0 이상
  # table 데이터의 variables 개수 입력
  colnames()
VNAM <- table %>% 
  filter(Group %in%"VNAM") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()

#####################################################################################

venn_df = list(`Control` = Control,
               `VNAM` = VNAM
)

#####################################################################################

#if (!require(devtools)) install.packages("devtools") 
#devtools::install_github("gaospecial/ggVennDiagram")

library(ggVennDiagram)

p <- ggVennDiagram(venn_df, 
                   label_alpha = 0,
                   category.names = c("Control", "VNAM"))+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  # ← 좌우 여백 늘리기
  scale_x_continuous(expand = expansion(mult = 0.1))


p
ggsave("venn4_Group.png", width = 7, height = 7, dpi = 500, 
       device = png, type = "cairo",  bg = "white")

#####################################################################################

#install.packages("UpSetR")
library(UpSetR)

table2 <- table[, c(11, 28:5416)]  %>%           # Group 열과 feature table만 추출 
  dplyr::group_by(Group) %>%               # Group 따라 feature값 합하기 
  dplyr::summarize_all(sum) %>% 
  tibble::column_to_rownames("Group")  %>% # rownames을 Group으로 설정
  t() %>%                                      # 전치행렬
  data.frame()

table2

#####################################################################################

# https://github.com/waldronlab/MicrobiomeWorkshop/blob/master/vignettes/MicrobiomeWorkshop.Rmd
table3 <- (table2>0) *1

p <- upset(data.frame(table3), 
           sets=colnames(data.frame(table3)),
           sets.bar.color = "#56B4E9",
           order.by = "freq",
           empty.intersections = "on")

p

png("venn5_Group.png", width = 10, height = 7, units = "in", res = 500,  type = "cairo",  bg = "white")
p
dev.off()

#####################################################################################
#####################################################################################

library(phyloseq)
library(dplyr)

rm(list = ls())

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")

ps <- readRDS("physeq_GutLungAxis_freq-min-1000.rds")

#####################################################################################

meta <- data.frame(sample_data(ps))
otu <- data.frame(otu_table(ps))

colnames(otu) <- gsub("\\.", "-", colnames(otu))
common <- intersect(rownames(meta), colnames(otu))
meta <- meta[common, , drop = FALSE]
otu  <- otu[, common, drop = FALSE]

table <- merge(meta, t(otu), by = "row.names")

#####################################################################################

Week0 <- table %>%
  filter(SamplingWeek %in%  "Week0") %>%
    dplyr::select(which(colSums(.[28:5416]) > 0)) %>%  # feature table이 column 6~643에 위치, 이 중 전체 합이 0 이상
  # table 데이터의 variables 개수 입력
  colnames()
Week1 <- table %>% 
  filter(SamplingWeek %in%"Week1") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()
Week2 <- table %>% 
  filter(SamplingWeek %in%"Week2") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()

#####################################################################################

venn_df = list(`Week0` = Week0,
               `Week1` = Week1,
               `Week2` = Week2
)

#####################################################################################

#if (!require(devtools)) install.packages("devtools") 
#devtools::install_github("gaospecial/ggVennDiagram")

library(ggVennDiagram)

p <- ggVennDiagram(venn_df, 
                   label_alpha = 0,
                   category.names = c("Week0", "Week1", "Week2"))+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  # ← 좌우 여백 늘리기
  scale_x_continuous(expand = expansion(mult = 0.1))


p
ggsave("venn4_SamplingWeek.png", width = 7, height = 7, dpi = 500, 
       device = png, type = "cairo",  bg = "white")

#####################################################################################

#install.packages("UpSetR")
library(UpSetR)

table2 <- table[, c(15, 28:5416)]  %>%           # Group 열과 feature table만 추출 
  dplyr::group_by(SamplingWeek) %>%               # Group 따라 feature값 합하기 
  dplyr::summarize_all(sum) %>% 
  tibble::column_to_rownames("SamplingWeek")  %>% # rownames을 Group으로 설정
  t() %>%                                      # 전치행렬
  data.frame()

table2

#####################################################################################

# https://github.com/waldronlab/MicrobiomeWorkshop/blob/master/vignettes/MicrobiomeWorkshop.Rmd
table3 <- (table2>0) *1

p <- upset(data.frame(table3), 
           sets=colnames(data.frame(table3)),
           sets.bar.color = "#56B4E9",
           order.by = "freq",
           empty.intersections = "on")

p

png("venn5_SamplingWeek.png", width = 10, height = 7, units = "in", res = 500,  type = "cairo",  bg = "white")
p
dev.off()

#####################################################################################
#####################################################################################

library(phyloseq)
library(dplyr)

rm(list = ls())

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")

ps <- readRDS("physeq_GutLungAxis_freq-min-1000.rds")

#####################################################################################

meta <- data.frame(sample_data(ps))
otu <- data.frame(otu_table(ps))

colnames(otu) <- gsub("\\.", "-", colnames(otu))
common <- intersect(rownames(meta), colnames(otu))
meta <- meta[common, , drop = FALSE]
otu  <- otu[, common, drop = FALSE]

table <- merge(meta, t(otu), by = "row.names")

#####################################################################################

Con_Week0 <- table %>%
  filter(Group.SamplingWeek %in%  "Con_Week0") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%  # feature table이 column 6~643에 위치, 이 중 전체 합이 0 이상
  # table 데이터의 variables 개수 입력
  colnames()
Con_Week1 <- table %>%
  filter(Group.SamplingWeek %in%  "Con_Week1") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%  # feature table이 column 6~643에 위치, 이 중 전체 합이 0 이상
  # table 데이터의 variables 개수 입력
  colnames()
Con_Week2 <- table %>%
  filter(Group.SamplingWeek %in%  "Con_Week2") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%  # feature table이 column 6~643에 위치, 이 중 전체 합이 0 이상
  # table 데이터의 variables 개수 입력
  colnames()

VNAM_Week0 <- table %>% 
  filter(Group.SamplingWeek %in%"VNAM_Week0") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()
VNAM_Week1 <- table %>% 
  filter(Group.SamplingWeek %in%"VNAM_Week1") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()
VNAM_Week2 <- table %>% 
  filter(Group.SamplingWeek %in%"VNAM_Week2") %>%
  dplyr::select(which(colSums(.[28:5416]) > 0)) %>%
  colnames()

#####################################################################################

venn_df = list(`Con_Week0` = Con_Week0,
               `Con_Week1` = Con_Week1,
               `Con_Week2` = Con_Week2,
               `VNAM_Week0` = VNAM_Week0,
               `VNAM_Week1` = VNAM_Week1,
               `VNAM_Week2` = VNAM_Week2
)

#####################################################################################

#if (!require(devtools)) install.packages("devtools") 
#devtools::install_github("gaospecial/ggVennDiagram")

library(ggVennDiagram)

p <- ggVennDiagram(venn_df, 
                   label_alpha = 0,
                   category.names = c("Con_Week0", "Con_Week1", "Con_Week2",
                                      "VNAM_Week0", "VNAM_Week1", "VNAM_Week2"))+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  # ← 좌우 여백 늘리기
  scale_x_continuous(expand = expansion(mult = 0.1))


p
ggsave("venn4_Group+SamplingWeek.png", width = 7, height = 7, dpi = 500, 
       device = png, type = "cairo",  bg = "white")

#####################################################################################

#install.packages("UpSetR")
library(UpSetR)

table2 <- table[, c(20, 28:5416)]  %>%           # Group 열과 feature table만 추출 
  dplyr::group_by(Group.SamplingWeek) %>%               # Group 따라 feature값 합하기 
  dplyr::summarize_all(sum) %>% 
  tibble::column_to_rownames("Group.SamplingWeek")  %>% # rownames을 Group으로 설정
  t() %>%                                      # 전치행렬
  data.frame()

table2

#####################################################################################

# https://github.com/waldronlab/MicrobiomeWorkshop/blob/master/vignettes/MicrobiomeWorkshop.Rmd
table3 <- (table2>0) *1

p <- upset(data.frame(table3), 
           sets=colnames(data.frame(table3)),
           sets.bar.color = "#56B4E9",
           order.by = "freq",
           empty.intersections = "on")

p

png("venn5_Group+SamplingWeek.png", width = 10, height = 7, units = "in", res = 500,  type = "cairo",  bg = "white")
p
dev.off()

#####################################################################################