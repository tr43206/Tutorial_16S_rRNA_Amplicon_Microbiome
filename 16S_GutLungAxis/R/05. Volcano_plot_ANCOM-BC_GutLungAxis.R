rm(list = ls())

library(readxl)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(mutoss)

setwd("C:/Users/tr432/Downloads/WeeksMerged_16sSeq/DA_testing/ANCOM-BC2")

df <- read_excel('ANCOMBC2_Group_primary_BKY_pseudo1_merged_heatmap용.xlsx', sheet = 'Week2')
head(df)

#####################################################################################

## (선택) 만약 BKY를 미리 보정 안했다면 -> BKY method FDR

## p_val을 숫자로(혹시 factor/char이면)
df$p_val <- as.numeric(df$p_mann)

## NA 제외하고 BKY 보정
keep <- !is.na(df$p_val)

bky_res <- mutoss::two.stage(pValues = df$p_val[keep],
                             alpha   = 0.05       # 필요시 다른 alpha로 변경
                             )

## q_val 열에 BKY q-value 저장
df$q_val <- NA_real_
df$q_val[keep] <- bky_res$adjPValues

#####################################################################################

## (선택) 만약 Genus만 미리 분리하지 않았다면

## 1) taxonomy 가 들어있는 열 이름 지정
tax_col <- "id"   # <-- 실제 계/문/강/목/과/속/종 정보가 들어있는 열 이름

## 2) ; 로 split해서 g__ 가 들어있는 부분(Genus)만 추출
df <- df %>%
  mutate(
    Genus_raw = sapply(strsplit(.data[[tax_col]], ";"), function(x) {
      x <- trimws(x)
      g_part <- x[grepl("g__", x)]
      if (length(g_part) == 0) return(NA_character_)
      g_part[1]
    }),
    ## 3) 앞의 "g__" 떼고 양쪽 공백 제거
    Genus = str_trim(str_remove(Genus_raw, "^g__"))
  ) %>%
  ## 4) Genus 가 비어있거나 NA 인 행 제거
  filter(!is.na(Genus), Genus != "")

## 이제 df에는 Genus 열이 생겼고,
## genus 정보가 없는 행은 제거된 상태

#####################################################################################

df$diffexpressed <- "NO"
df$diffexpressed[df$lfc > 1 & df$q_val < 0.05] <- "UP"  #Edit: lfc, p_val(or q-val) 조정
df$diffexpressed[df$lfc < -1 & df$q_val < 0.05] <- "DOWN"  #Edit: lfc, p_val(or q-val) 조정
head(df[order(df$q_val) & df$diffexpressed == 'DOWN', ])

df_qfilter <- df[df$q_val < 0.05,]

# Genus를 문자로 고정(공백도 제거)
df$Genus <- trimws(as.character(df$Genus))


## UP/DOWN 전체에서 Top 20 뽑는거
#df$delabel <- ifelse(df$id %in% head(df[order(df$q_val), "id"], 20), df$id, NA)  #Edit: order, 숫자(top 몇개까지 annotation할지) 조정


## (선택1) 각 그룹에서 q-value 기준 상위 n개씩만 annotation 할때
up_top  <- head(df[df$diffexpressed == "UP", ][order(df$lfc[df$diffexpressed == "UP"], decreasing = TRUE),  "Genus"], 10)  #Edit
down_top <- head(df[df$diffexpressed == "DOWN", ][order(df$lfc[df$diffexpressed == "DOWN"]), "Genus"], 30)  #Edit


## (선택2) 각 그룹에서 q-value < 0.05인거 전부 annotation 할때
up_idx   <- which(df$diffexpressed == "UP")
down_idx <- which(df$diffexpressed == "DOWN")

up_top   <- df$Genus[up_idx][order(df$q_val[up_idx])]
down_top <- df$Genus[down_idx][order(df$q_val[down_idx])]

ids_to_label <- unique(c(up_top, down_top))

df$delabel <- ifelse(df$Genus %in% ids_to_label, df$Genus, NA_character_)

sum(!is.na(df$delabel))   # 0 아니면 라벨 준비 완료

#ids_to_label <- c(up_top, down_top)
#df$delabel <- ifelse(df$Genus %in% ids_to_label, df$Genus, NA)  #Edit
#head(df[order(df$q_val), "Genus"], 20)  #Edit


## (선택3) 특정 균주만 annotation
df$delabel <- ifelse(df$Genus == 'Phocaeicola_A', df$Genus, NA)  #Edit


theme_set(theme_classic(base_size = 13) +  #Edit: base_size 조정
            theme(
              axis.title.y = element_text(face = "bold", margin = ggplot2::margin(0,20,0,0), size = rel(1.1), color = 'black'),  #Edit: rel 조정
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(20,0,0,0), size = rel(1.1), color = 'black'),  #Edit: rel 조정
              plot.title = element_text(hjust = 0.5)  #Edit: hjust 조정
            ))

volcano <- ggplot(data = df, aes(x = lfc, y = -log10(q_val), col = diffexpressed, label = delabel)) +  #Edit: -log10() p_val로 할지 q_val로 할지 결정해야됨
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = 'dashed') +  #Edit: lfc 기준 바꿀때 변경
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("firebrick", "grey", "royalblue"),
                     labels = c("Control", "Not significant", "VNAM")) +
  coord_cartesian(ylim = c(0, 7), xlim = c(-6, 6)) +  #Edit: ylim, xlim 조정
  labs(color = 'Group',
       x = expression("log"[2]*"FC (VNAM/Control)"), y = expression("-log"[10]*"(q-value)")) +  #Edit: -log p-value로 할지 q-value로 할지 결정해야됨
  scale_x_continuous(breaks = seq(-6, 6, 1)) +  #Edit: seq 조정(최소값, 최대값, 간격)
  scale_y_continuous(breaks = seq(0, 7, 1)) +  #Edit: seq 조정(최소값, 최대값, 간격)
  ggtitle('ANCOM-BC2 (Week2)') +  #Edit: ggtitle 수정
  geom_text_repel(max.overlaps = Inf) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0))  ##**##

print(volcano)

ggsave("C:/Users/tr432/Downloads/filter_GutLungAxis_Week2/ancombc2_Group_pseudo1_volcano.png", width = 10, height = 10, dpi = 500)
#write.csv(df, file = "filter_GutLungAxis_Week2/ancombc2_BKY.csv", row.names = FALSE)

#pdf(file = "group_volcano_g.pdf", width = 8, height = 6)
#volcano
#dev.off()

#####################################################################################

