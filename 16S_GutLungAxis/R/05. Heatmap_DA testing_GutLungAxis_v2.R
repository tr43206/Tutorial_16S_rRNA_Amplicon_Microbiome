rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd('C:/Users/tr432/Downloads/WeeksMerged_16sSeq/DA_testing/MaAsLin3')

################################################################################
xlsx_path <- 'maaslin3_Group+Initial_b.w_prev_BKY_merged_heatmap용.xlsx'
alpha_q   <- 0.05   # 유의성 기준(q-val)
digits_eff <- 2     # 타일에 표시할 effect size 소수점

################################################################################
# =========================
# 1) 시트 읽기 + week 라벨 부여
# =========================
read_week <- function(sheet, week_label){
  read_excel(xlsx_path, sheet = sheet) %>%
    transmute(
      week  = week_label,
      taxon = as.character(Genus),
      eff   = as.numeric(coef),  #Edit: lda or effect or lfc or coef
      p_val = as.numeric(p_val),
      q_val = as.numeric(q_val)
    )
}

df0 <- read_week('Week0', '0w')
df1 <- read_week('Week1', '1w')
df2 <- read_week('Week2', '2w')

df_all <- bind_rows(df0,
                    df1, df2)

################################################################################
# =========================
# 2) y축(Genus) 선정: 어느 주차든 q < alpha 인 Genus '합집합'
# =========================
tax_keep <- df_all %>%
  filter(!is.na(q_val), q_val < alpha_q) %>%
  distinct(taxon) %>%
  pull(taxon)

# (옵션) 만약 "3개 주차 모두에서 유의"한 교집합을 원하면 아래로 교체
# tax_keep <- df_all %>%
#   filter(!is.na(q_val), q_val < alpha_q) %>%
#   count(taxon) %>%
#   filter(n == 3) %>%
#   pull(taxon)

# keep이 비면(유의한 genus 없음) 전체를 쓰거나 중단
if (length(tax_keep) == 0) stop('q_val < alpha 기준을 만족하는 Genus가 없습니다.')

################################################################################
# =========================
# 3) Plot용 long 데이터: 유의하면 eff, 아니면 0 / 텍스트 색상도 같이
# =========================
sig_color <- '#006400'

df_fig <- df_all %>%
  filter(taxon %in% tax_keep) %>%
  mutate(
    week = factor(week, levels = c('0w', '1w', '2w')),
    sig  = !is.na(q_val) & (q_val < alpha_q),
#    value = ifelse(sig, round(eff, digits_eff), 0),
    value = round(eff, digits_eff),
    
    label = formatC(value, format = "f", digits = digits_eff),
    
    star = case_when(
      !is.na(q_val) & q_val < 0.0001 ~ '****',
      !is.na(q_val) & q_val < 0.001  ~ '***',
      !is.na(q_val) & q_val < 0.01   ~ '**',
      !is.na(q_val) & q_val < 0.05   ~ '*',
      TRUE ~ ''
    ),
    
    text_color = ifelse(sig, sig_color, 'black'),
    fontface = ifelse(sig, "bold", "plain")
  ) %>%
  dplyr::select(taxon, week, value, q_val, sig, label, star, text_color, fontface)

# y축 정렬(기본: 이름순). 원하면 효과크기 기준으로 정렬 가능.
df_fig <- df_fig %>%
  mutate(taxon = factor(taxon, levels = sort(unique(taxon))))

################################################################################
# =========================
# 4) 색상 범위 설정 (midpoint=0)
# =========================
lo  <- floor(min(df_fig$value, na.rm = TRUE))
up  <- ceiling(max(df_fig$value, na.rm = TRUE))
mid <- 0

################################################################################
# =========================
# 5) Heatmap
# =========================
n_x <- length(levels(df_fig$week))
n_y <- length(levels(df_fig$taxon))
line_w <- 0.35

p <- ggplot(df_fig, aes(x = week, y = taxon, fill = value)) +
  geom_tile(color = NA) +   # 타일 자체 테두리는 제거
  
  # 안쪽 + 바깥쪽 테두리를 동일 두께로 직접 그림
  geom_vline(
    xintercept = seq(0.5, n_x + 0.5, by = 1),
    color = 'black',
    linewidth = line_w
  ) +
  geom_hline(
    yintercept = seq(0.5, n_y + 0.5, by = 1),
    color = 'black',
    linewidth = line_w
  ) +
  
  scale_fill_gradient2(
    low = '#B22222', mid = "white", high = '#4169E1',
    midpoint = mid, limits = c(lo, up),
    na.value = 'white', name = 'β coefficient'  #Edit: 'LDA score' or 'Effect size' or 'Log2FC' or 'β coefficient'
  ) +
  
  geom_text(
    aes(label = label, color = text_color, fontface = fontface),
    size = 4
  ) +
  
  geom_text(
    data = df_fig %>% filter(sig),
    aes(label = star),
    color = sig_color,
    fontface = "bold",
    size = 6.5,
    nudge_x = 0.25,
    nudge_y = 0.09  #Edit: 0.09 단위로 조절하는걸 추천
  ) +
  
  scale_color_identity(guide = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  
  labs(x = NULL, y = NULL, title = 'MaAsLin3 (prevalence)') +  #Edit: 'LEfSe' or 'ALDEx2' or 'ANCOM-BC2' or 'MaAsLin3 (abundance)' or 'MaAsLin3 (prevalence)
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    # 중요: panel.border 제거
    panel.border = element_blank(),
    
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.6, "cm")
  )

p


ggsave('maaslin3_prev_Group+Initial_b.w_uncultured 제외_v2.png', p, width = 9, height = 10, dpi = 500)
