# I compare the UKBB trajectory results with the Danish results
library(dplyr)
library(stringr)
library(MatchIt)
library(plyr)
library(tibble)
library(igraph)
library(survival)
library(Matrix)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(reshape2)
library(viridis)
library(readxl)
library(stringr)
library(data.table)
library(glmnet)
library(reticulate)
library(tidyverse)
library(foreach)
library(doParallel)
rm(list=ls())

# output folder
out.dir = "output/ICD10_level3_diseases/trajectory_vali"
dir.create(out.dir, recursive = T,showWarnings = FALSE)
input.dir = 'output/ICD10_level3_diseases'
comm_file.dir = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB'

location = 'mac'
if (location == 'tetralith') {
  setwd("/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3")
  comm_file.dir = '/proj/spatial_pre_mal/users/x_yelzh/Projects/Common files/UKBB' # tetralith
}

################## 
#1 Load UKBB results ####
final_results_df <- readRDS(file = paste0(out.dir,'/Step-RR_direction_results_doi.RDS'))
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'
# scp -r '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali' x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'  

final_trajectories <- final_results_df %>% mutate(
  p.value.RR.adj = p.adjust(p_val_RR, method = "fdr"),
  RR_significant = ifelse(p.value.RR.adj < 0.05, '1', '-1')) %>%
  filter(RR_significant == '1') %>% 
  mutate(p_val_direction.adj = p.adjust(p_val_direction, method = "fdr"),
         direction_yes_no = ifelse(p_val_direction.adj < 0.05, '1', '-1')) %>% 
  filter(direction_yes_no == '1') %>%
  arrange(D2,desc(RR))

colnames(final_trajectories)

#2 Load Danish results ####
danish_study_diseases <- read.csv2(paste0('input/DanishReg_NatComm2020_PMID33009368/RR_from_DanishReg_NatComm2020_PMID33009368.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
danish_study_diseases <- danish_study_diseases %>% filter(D2 %in% final_trajectories$D2) %>%
  filter(RR > 1) %>% filter(direction_yes_no == '1') 
colnames(danish_study_diseases)

#3 Compare UKBB and Danish results, keep all unique pairs####
merged_trajectory_results <- final_trajectories %>%
  select(D1, D2, D1_name,D2_name, RR, direction_yes_no) %>% 
  inner_join(danish_study_diseases %>% select(D1, D2, RR_Danish = RR, direction_yes_no_Danish = direction_yes_no),
            by = c("D1", "D2")) %>% mutate(RR = as.numeric(RR),
                                             RR_Danish = as.numeric(RR_Danish), D1_D2_pairs = paste0(D1,'-',D2))

#4 Plot comparison of RR ####
library(ggrepel)
cor_test <- cor.test(log2(merged_trajectory_results$RR_Danish), log2(merged_trajectory_results$RR), method = "spearman")
cor_label <- sprintf("Spearman R: %.3f\np-value: %.3e", 
                     cor_test$estimate, 
                     cor_test$p.value)
p = ggplot(merged_trajectory_results, aes(x = log2(RR_Danish), y = log2(RR), label = D1_D2_pairs)) +
  geom_point(color = 'blue', alpha = 0.5, size = 2) +
  geom_smooth(method = 'lm', color = 'grey', se = FALSE) +
  labs(title = 'Comparison of RR between UKBB and Danish Study',
       x = 'Log2 RR (Danish Study)',
       y = 'Log2 RR (UKBB)') + 
  geom_text_repel(size = 3, max.overlaps = 20) +
  annotate("text", 
           x = min(log2(merged_trajectory_results$RR_Danish), na.rm = TRUE), 
           y = max(log2(merged_trajectory_results$RR), na.rm = TRUE), 
           label = cor_label, 
           hjust = 0, vjust = 1, # Aligns text to the top-left corner
           size = 4, fontface = "bold", color = "black") +
  theme_classic()  
ggsave(p, file = paste0(out.dir,'/RR_comparison_UKBB_Danish_study.pdf'), width = 6, height = 5)

# ## Plot for each D2 disease ####
# unique_D2s <- unique(merged_trajectory_results$D2_name)
# for (d2 in unique_D2s) {
#   subset_data <- merged_trajectory_results %>% filter(D2_name == d2)
#   
#   cor_test_subset <- cor.test(log2(subset_data$RR_Danish), log2(subset_data$RR), method = "spearman")
#   cor_label_subset <- sprintf("Spearman R: %.3f\np-value: %.3e", 
#                               cor_test_subset$estimate, 
#                               cor_test_subset$p.value)
#   
#   p_subset = ggplot(subset_data, aes(x = log2(RR_Danish), y = log2(RR), label = D1_name)) +
#     geom_point(color = 'blue', alpha = 0.5, size = 2) +
#     geom_smooth(method = 'lm', color = 'grey', se = FALSE) +
#     labs(title = paste0('Comparison of RR for D2: ', d2),
#          x = 'Log2 RR (Danish Study)',
#          y = 'Log2 RR (UKBB)') + 
#     geom_text_repel(size = 3, max.overlaps = 20) +
#     annotate("text", 
#              x = min(log2(subset_data$RR_Danish), na.rm = TRUE), 
#              y = max(log2(subset_data$RR), na.rm = TRUE), 
#              label = cor_label_subset, 
#              hjust = 0, vjust = 1, # Aligns text to the top-left corner
#              size = 4, fontface = "bold", color = "black") +
#     theme_classic()  
#   
#   ggsave(p_subset, file = paste0(out.dir,'/RR_comparison_UKBB_Danish_study_D2_', d2, '.pdf'), width = 6, height = 5)
# }
