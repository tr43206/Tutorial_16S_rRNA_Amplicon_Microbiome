sessionInfo()

#####################################################################################
rm(list = ls())

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("biobakery/maaslin3")
library(maaslin3)
library(mutoss)

#####################################################################################
### 3. Microbiome association detection with MaAsLin 3

## 3.1 MaAsLin 3 input

setwd("C:/Users/tr432/Downloads/Figures_GutLungAxis/yes0")

### 상대적 풍부도 데이터
# Read abundance table
taxa_table_name <- "rel_gg2_genus_table.tsv"
taxa_table <- read.csv(taxa_table_name, sep = '\t', row.names = 1,
                       skip = 1, check.names = FALSE, comment.char = "")
taxa_table <- t(taxa_table)

# Read metadata table
metadata_name <- "RUN45_Metadata_GutLungAxis.txt"
metadata <- read.csv(metadata_name, sep = '\t', row.names = 1)


# Factor the categorical variables to test IBD against healthy controls
metadata$SamplingWeek <-
  factor(metadata$SamplingWeek, levels = c('Week0', 'Week1', 'Week2'))
metadata$Group <-
  factor(metadata$Group, levels = c('Control', 'VNAM'))

colnames(taxa_table) <- gsub("\\.", "-", colnames(taxa_table))

taxa_table[1:5, 1:5]
metadata[1:5, 1:5]

#####################################################################################

set.seed(1)
fit_out <- maaslin3(
  input_data    = taxa_table,
  input_metadata = metadata,
  output        = "maaslin3_SamplingWeek+Group+BM_wbc+b.w_BH_rel",  #Edit
  formula       = "~ SamplingWeek + Group + BM_wbc + b.w",  #Edit
  normalization = "TSS",
  transform     = "LOG",
  correction    = "BH",      # 혹은 "none" 후 BKY 수동 보정
  augment       = TRUE,
  standardize   = TRUE,
  max_significance = 0.1,
  median_comparison_abundance  = TRUE,  ## 상대적 풍부도엔 써도되지만 절대 풍부도엔 쓰면 안됨
  median_comparison_prevalence = FALSE,
#  max_pngs     = 20,
  cores        = 12
)


#scatter_plots <- maaslin_plot_results_from_output(
#  output = "../GutLungAxis_maaslin3_Group_BH_rel",
#  metadata = metadata,
#  normalization = "TSS",
#  transform = "LOG",
#  median_comparison_abundance  = TRUE,
#  median_comparison_prevalence = FALSE,
#  max_significance = 0.05,
#  max_pngs = 20
#)


res_abund <- fit_out$fit_data_abundance$results %>%
  dplyr::select(-qval_individual) %>%
  dplyr::select(-qval_joint)

keep_ind <- !is.na(res_abund$pval_individual)
bky_ind  <- mutoss::two.stage(res_abund$pval_individual[keep_ind], alpha = 0.05)
res_abund$q_ind_bky <- NA_real_
res_abund$q_ind_bky[keep_ind] <- bky_ind$adjPValues

keep_joint <- !is.na(res_abund$pval_joint)
bky_joint  <- mutoss::two.stage(res_abund$pval_joint[keep_joint], alpha = 0.05)
res_abund$q_joint_bky <- NA_real_
res_abund$q_joint_bky[keep_joint] <- bky_joint$adjPValues

write.csv(res_abund, "maaslin3_SamplingWeek+Group+BM_wbc+b.w_BKY_rel.csv", row.names = FALSE)  #Edit

#####################################################################################
## 3.2 Running MaAsLin 3

set.seed(1)
fit_out <- maaslin3(input_data = taxa_table,
                    input_metadata = metadata,
                    output = 'GutLungAxis_maaslin3_BKY',
                    formula = '~ Group',
                    normalization = 'TSS',
                    transform = 'LOG',
                    correction = "BH",
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 20,
                    cores = 12)


res_abund <- fit_out$fit_data_abundance$results

## individual p-values에 대한 BKY
keep_ind <- !is.na(res_abund$pval_individual)
bky_ind <- mutoss::two.stage(pValues = res_abund$pval_individual[keep_ind],
                             alpha = 0.05)
res_abund$q_ind_bky <- NA_real_
res_abund$q_ind_bky[keep_ind] <- bky_ind$adjPValues

## joint p-values에 대한 BKY
keep_joint <- !is.na(res_abund$pval_joint)
bky_joint <- mutoss::two.stage(pValues = res_abund$pval_joint[keep_joint],
                               alpha = 0.05)
res_abund$q_joint_bky <- NA_real_
res_abund$q_joint_bky[keep_joint] <- bky_joint$adjPValues

write.csv(res_abund, "GutLungAxis_maaslin3_BKY.csv", row.names = FALSE)

#####################################################################################
## 3.3 MaAsLin 3 output

# This section is necessary for updating the
# summary plot and the association plots

### 그룹명, 개체명 바꾸는 단계
# Rename results file with clean titles
all_results <- read.csv('hmp2_output/all_results.tsv', sep='\t')
all_results <- all_results %>%
  mutate(metadata = case_when(metadata == 'age' ~ 'Age',
                              metadata == 'antibiotics' ~ 'Abx',
                              metadata == 'diagnosis' ~ 'Diagnosis',
                              metadata == 'dysbiosis_state' ~ 'Dysbiosis',
                              metadata == 'reads' ~ 'Read depth'),
         value = case_when(value == 'dysbiosis_CD' ~ 'CD',
                           value == 'dysbiosis_UC' ~ 'UC',
                           value == 'Yes' ~ 'Used', # Antibiotics
                           value == 'age' ~ 'Age',
                           value == 'reads' ~ 'Read depth',
                           TRUE ~ value),
         feature = gsub('_', ' ', feature) %>%
           gsub(pattern = 'sp ', replacement = 'sp. '))

# Write results
write.table(all_results, 'hmp2_output/all_results.tsv', sep='\t')

# Set the new heatmap and coefficient plot variables and order them
heatmap_vars = c('Dysbiosis UC', 'Diagnosis UC',
                 'Abx Used', 'Age', 'Read depth')
coef_plot_vars = c('Dysbiosis CD', 'Diagnosis CD')

# This section is necessary for updating the association plots
taxa_table_copy <- taxa_table
colnames(taxa_table_copy) <- gsub('_', ' ', colnames(taxa_table_copy)) %>%
  gsub(pattern = 'sp ', replacement = 'sp. ')

# Rename the features in the norm transformed data file
data_transformed <-
  read.csv('hmp2_output/features/data_transformed.tsv', sep='\t')
colnames(data_transformed) <-
  gsub('_', ' ', colnames(data_transformed)) %>%
  gsub(pattern = 'sp ', replacement = 'sp. ')
write.table(data_transformed,
            'hmp2_output/features/data_transformed.tsv',
            sep='\t', row.names = FALSE)

# Rename the metadata like in the outputs table
metadata_copy <- metadata
colnames(metadata_copy) <-
  case_when(colnames(metadata_copy) == 'age' ~ 'Age',
            colnames(metadata_copy) == 'antibiotics' ~ 'Abx',
            colnames(metadata_copy) == 'diagnosis' ~ 'Diagnosis',
            colnames(metadata_copy) == 'dysbiosis_state' ~ 'Dysbiosis',
            colnames(metadata_copy) == 'reads' ~ 'Read depth',
            TRUE ~ colnames(metadata_copy))
metadata_copy <- metadata_copy %>%
  mutate(Dysbiosis = case_when(Dysbiosis == 'dysbiosis_UC' ~ 'UC',
                               Dysbiosis == 'dysbiosis_CD' ~ 'CD',
                               Dysbiosis == 'none' ~ 'None') %>%
           factor(levels = c('None', 'UC', 'CD')),
         Abx = case_when(Abx == 'Yes' ~ 'Used',
                         Abx == 'No' ~ 'Not used') %>%
           factor(levels = c('Not used', 'Used')),
         Diagnosis = case_when(Diagnosis == 'nonIBD' ~ 'non-IBD',
                               TRUE ~ Diagnosis) %>%
           factor(levels = c('non-IBD', 'UC', 'CD')))

# Recreate the plots
scatter_plots <- maaslin_plot_results_from_output(
  output = 'hmp2_output',
  metadata = metadata_copy,
  normalization = "TSS",
  transform = "LOG",
  median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,
  max_significance = 0.1,
  max_pngs = 20)

#####################################################################################

