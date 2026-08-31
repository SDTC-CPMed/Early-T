# In this version, I used the same method for Danish data mentioned in https://doi.org/10.1038/s41467-020-18682-4
# But instead of using their results, I run the filtering pipeline again to get RR and directionality p-values in UKBB
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
#1 Load diseases in UKBB ####
disease_n_patient <- readRDS(file = paste0(input.dir,'/disease_n_patient.RDS'))
disease_n_patient_Prot <- readRDS(file = paste0(input.dir,'/disease_n_patient_proteomics.RDS'))

### Filter1. Only keep diseases with at least n_patient_threshold number of patients ####
n_patient_threshold = 100
n_patient_prot_threshold = 80
disease_n_patient_filtered = disease_n_patient %>% filter(Number_pt >= n_patient_threshold)
disease_n_patient_Prot_filtered = disease_n_patient_Prot %>% filter(Number_pt >= n_patient_prot_threshold)

All_dis_pairs <- expand.grid(D1 = disease_n_patient_filtered$ICD10,
                             D2 = disease_n_patient_filtered$ICD10,
                             stringsAsFactors = FALSE) %>%  filter(D1 != D2) #1,438,800 pairs

#2 Load diseases from the Danish study (DanishReg_NatComm2020_PMID33009368) ####
danish_study_diseases <- read.csv2(paste0('input/DanishReg_NatComm2020_PMID33009368/RR_from_DanishReg_NatComm2020_PMID33009368.tsv'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
head(danish_study_diseases)

# add disease name for D1 and D2 using icd10_level3_disease_names
danish_study_diseases <- danish_study_diseases %>%
  left_join(icd10_level3_disease_names, by = c("D1" = "coding")) %>%
  mutate(D1_name = meaning) %>% select(-meaning) %>%
  left_join(icd10_level3_disease_names, by = c("D2" = "coding")) %>%
  mutate(D2_name = meaning) %>% select(-meaning) %>%
  filter(direction_yes_no == 1)
head(danish_study_diseases)
write_xlsx(danish_study_diseases, path = paste0('input/DanishReg_NatComm2020_PMID33009368/RR_from_DanishReg_NatComm2020_PMID33009368_sig_ICDnames.xlsx'))

# filter p.value.direction < 0.05 and shared with UKBB and to get all unique D1 and D2 disease pairs
danish_study_diseases_pairs_filtered <- danish_study_diseases %>%
  filter(direction_yes_no == 1) %>% 
  filter((D1 %in% disease_n_patient_filtered$ICD10 ) & (D2 %in% disease_n_patient_filtered$ICD10)) %>% mutate(D1_D2_pairs = paste0(D1, "_", D2)) %>% unique
head(danish_study_diseases_pairs_filtered)
dim(danish_study_diseases_pairs_filtered)

# filter p.value.direction < 0.05 and pull all unique D1 and D2 diseases
danish_study_diseases_filtered <- danish_study_diseases_pairs_filtered %>%
  select(D1, D2) %>%
  pivot_longer(cols = c(D1, D2), names_to = "D_type", values_to = "ICD10_Disease") %>%
  select(ICD10_Disease) %>% pull(ICD10_Disease) %>%
  unique() 
length(danish_study_diseases_filtered)

#3 Load UKBB clinic data ####
All_clinic = readRDS(paste0(comm_file.dir,'/ukb_clinic.rda')) 
All_clinic = All_clinic %>% 
  mutate(across(c(date_of_attending_assessment_centre_f53_0_0, date_of_death_f40000_0_0, date_of_death_f40000_1_0, date_lost_to_followup_f191_0_0), as.Date)) %>%
  mutate(assessment_year = as.numeric(format(date_of_attending_assessment_centre_f53_0_0, "%Y"))) %>%
  mutate(assessment_month = as.numeric(format(date_of_attending_assessment_centre_f53_0_0, "%m"))) %>%
  mutate(sex_f31_0_0 = as.factor(sex_f31_0_0),
         max_age = case_when(age_at_recruitment_f21022_0_0 < 20 ~ '20-',
           age_at_recruitment_f21022_0_0 >= 20 & age_at_recruitment_f21022_0_0 < 30 ~ '20-29',
           age_at_recruitment_f21022_0_0 >= 30 & age_at_recruitment_f21022_0_0 < 40 ~ '30-39',
           age_at_recruitment_f21022_0_0 >= 40 & age_at_recruitment_f21022_0_0 < 50 ~ '40-49',
           age_at_recruitment_f21022_0_0 >= 50 & age_at_recruitment_f21022_0_0 < 60 ~ '50-59',
           age_at_recruitment_f21022_0_0 >= 60 & age_at_recruitment_f21022_0_0 < 70 ~ '60-69',
           age_at_recruitment_f21022_0_0 >= 70 ~ '70+',
           TRUE ~ NA_character_
         ),
         max_age = as.factor(max_age, levels = c('30-39', '40-49', '50-59', '60-69', '70+')),
         ethnic_background_f21000_0_0 = as.factor(ethnic_background_f21000_0_0),
         delta_diag_enroll = NA,
         OSstatus.raw = ifelse(is.na(date_of_death_f40000_0_0) & is.na(date_of_death_f40000_1_0), 0, 1),
         visible = ifelse(OSstatus.raw == 1, 
                          (date_of_death_f40000_0_0 - date_of_attending_assessment_centre_f53_0_0)/365,
                          ifelse(is.na(date_lost_to_followup_f191_0_0),
                                 (as.Date("2022-12-30") - date_of_attending_assessment_centre_f53_0_0)/365,
                                 (pmin(date_lost_to_followup_f191_0_0, "2022-12-30") - date_of_attending_assessment_centre_f53_0_0)/365))) %>% 
  filter(!is.na(ethnic_background_f21000_0_0),
         !is.na(age_at_recruitment_f21022_0_0),
         !is.na(uk_biobank_assessment_centre_f54_0_0),
         !is.na(sex_f31_0_0)) %>% 
  mutate(eid = as.character(eid))

#4 Main analysis  ####
# 
# D2_icd = "C18"
# D2_name = disease_n_patient %>% filter(ICD10 == D2_icd) %>% pull(ICD10_Disease)
# message('Analyzing end disease: ', D2_name)
# # get the all D1 diseases that have trajectory to D2 in Danish study
# all_D1 <- danish_study_diseases_pairs_filtered %>% filter(D2 == D2_icd) %>% pull(D1) %>% unique()
# D2_data = readRDS(paste0(input.dir,'/ukb_Z33_Pregnant_state_incidental_subset_clinic.rda')) 
# 
# # danish_study_diseases_pairs_filtered = danish_study_diseases_pairs_filtered %>% filter(D2 == D2_icd)

library(dplyr)
library(MatchIt)
## Configuration ####
MIN_EXPOSED <- 100
PREFILTER_P_THRESHOLD <- 0.05 / 1438800 # Bonferroni corrected threshold (example: 0.05 / 1,438,800 potential pairs)

## A. Pre-filter ####  
prefilter_results <- list()
for(i in 1:nrow(All_dis_pairs)) { #nrow(All_dis_pairs)
  d1_code <- All_dis_pairs$D1[i]
  d2_code <- All_dis_pairs$D2[i]
  D1_name = disease_n_patient %>% filter(ICD10 == d1_code) %>% pull(ICD10_Disease)
  D2_name = disease_n_patient %>% filter(ICD10 == d2_code) %>% pull(ICD10_Disease)
  
  message('Processing disease pair: ',i, ':', D1_name, ' -> ', D2_name)
  
  ### Load D1 and D2 specific data  ####
  file_d1 <- list.files(input.dir,pattern = paste0("ukb_", d1_code,"_"), full.names = TRUE)
  file_d2 <- list.files(input.dir,pattern = paste0("ukb_", d2_code,"_"), full.names = TRUE)
  
  df_d1 <- readRDS(file_d1)  %>% mutate(eid = as.character(eid))
  df_d2 <- readRDS(file_d2)  %>% mutate(eid = as.character(eid))
  
  ### Filter2. Skip if less than 20 overlap eids between D1 and D2  ####
  if (length(intersect(df_d1$eid, df_d2$eid)) < 20) {
    prefilter_results[[i]] <- data.frame(
      D1 = d1_code,
      D2 = d2_code,
      prefilter.p = NA,
      D1_counts = nrow(df_d1),
      D2_counts = nrow(df_d2),
      D1_name = D1_name,
      D2_name = D2_name,Note = 'Less than 20 comorbidity cases')
    
    message("Less than 20 comorbidity cases, skipping")
    next
  }
  
  ### Filter3. Bernoulli Trial ####
  # Probability of D2 in the general population (Prevalence)
  p_d2_population <- nrow(df_d2) / nrow(All_clinic)
  
  # Join to find any sequence where D1 occurred before D2
  timing_df <- inner_join(
    df_d1 %>% select(eid, diagnosis_time_d1 = diagnosis_time), 
    df_d2 %>% select(eid, diagnosis_time_d2 = diagnosis_time), 
    by = "eid"
  )
  
  # Check if D1 happened before D2
  c_exposed_lifetime <- timing_df %>% 
    filter(diagnosis_time_d2 > diagnosis_time_d1) %>% 
    nrow()
  
  n_exposed <- nrow(df_d1)
  
  # Bernoulli trial: Is the lifetime progression rate higher than general prevalence?
  prefilter_test <- binom.test(c_exposed_lifetime, n_exposed, p = p_d2_population, alternative = "greater")
  
    prefilter_results[[i]] <- data.frame(
      D1 = d1_code,
      D2 = d2_code,
      prefilter.p = prefilter_test$p.value,
      D1_counts = nrow(df_d1),
      D2_counts = nrow(df_d2),
      D1_name = D1_name,
      D2_name = D2_name,Note = NA)
}
saveRDS(prefilter_results, file = paste0(out.dir,'/Step-Prefilter_intermediate_results.RDS'))
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'
# scp -r '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali' x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'  

prefilter_results <- readRDS(file = paste0(out.dir,'/Step-Prefilter_intermediate_results.RDS')) # 1438800
prefilter_df <- bind_rows(prefilter_results)
prefilter_df <- prefilter_df %>% mutate(prefilter.p.adjust = p.adjust(prefilter.p, method = "bonferroni"))

saveRDS(prefilter_df, file = paste0(out.dir,'/Step-Prefilter_result.RDS'))
write.csv2(prefilter_df, file = paste0(out.dir,'/Step-Prefilter_result.csv'), row.names = FALSE)

prefilter_df <- readRDS(file = paste0(out.dir,'/Step-Prefilter_result.RDS'))



