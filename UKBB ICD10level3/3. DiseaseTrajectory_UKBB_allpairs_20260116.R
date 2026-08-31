# In this version, I used the same method for Danish data mentioned in https://doi.org/10.1038/s41467-020-18682-4
# But instead of using their results, I run the whole pipe line again to get RR and directionality p-values in UKBB
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
  
  # ### Filter1. D1 and D2 should have at least MIN_EXPOSED number of patients ####
  # if (nrow(df_d1) < MIN_EXPOSED & nrow(df_d2) < MIN_EXPOSED) {
  #   prefilter_results[[i]] <- data.frame(
  #     D1 = d1_code,
  #     D2 = d2_code,
  #     prefilter.p = NA,
  #     RR = NA,
  #     p.value.direction = NA,
  #     D1_counts = nrow(df_d1),
  #     D2_counts = nrow(df_d2),
  #     D1_name = D1_name,
  #     D2_name = D2_name,Note = 'Less than 100 cases in D1 or D2')
  #   
  #   message("Too few cases, skipping")
  #   next
  # }
  
  ### Filter2. Skip if less than 20 overlap eids between D1 and D2  ####
  if (length(intersect(df_d1$eid, df_d2$eid)) < 20) {
    prefilter_results[[i]] <- data.frame(
      D1 = d1_code,
      D2 = d2_code,
      prefilter.p = NA,
      RR = NA,
      p.value.direction = NA,
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
      RR = NA,
      p.value.direction = NA,
      D1_counts = nrow(df_d1),
      D2_counts = nrow(df_d2),
      D1_name = D1_name,
      D2_name = D2_name,Note = NA)
}
saveRDS(prefilter_results, file = paste0(out.dir,'/Step-Prefilter_intermediate_results.RDS'))
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'

prefilter_results <- readRDS(file = paste0(out.dir,'/Step-Prefilter_intermediate_results.RDS'))
prefilter_results[[1]]
prefilter_df <- bind_rows(prefilter_results)
prefilter_df <- prefilter_df %>% mutate(prefilter.p.adjust = p.adjust(prefilter.p, method = "bonferroni"))

saveRDS(prefilter_df, file = paste0(out.dir,'/Step-Prefilter_result.RDS'))
write.csv2(prefilter_df, file = paste0(out.dir,'/Step-Prefilter_result.csv'), row.names = FALSE)
prefilter_df <- readRDS(file = paste0(out.dir,'/Step-Prefilter_result.RDS'))


# subset trajectories to end diseases of interest ####
disease_oi <- c("C50", # Breast C
                "C61", # Prostate C
                "C18", # Colon C
                "C19", #Malignant_neoplasm_of_rectosigmoid_junction
                "C20",#Malignant_neoplasm_of_rectum"
                "C43", #C43_Malignant_melanoma_of_skin
                "C53",#C53_Malignant_neoplasm_of_cervix_uteri
                "C54",#Malignant_neoplasm_of_corpus_uteri
                "C55",#Malignant_neoplasm_of_uterus_part_unspecified
                "C62",#Malignant_neoplasm_of_testis
                "C56",#Malignant_neoplasm_of_ovary
                "C67",#Malignant_neoplasm_of_bladder
                "C64",#Malignant_neoplasm_of_kidney_except_renal_pelvis
                "C34",#Malignant_neoplasm_of_bronchus_and_lung
                "C25" #Malignant_neoplasm_of_pancreas
                ) #
sex_specific_diseases <- readRDS(file = paste0(input.dir,'/sex_specific_diseases.rds'))

# prefilter_df  <- readRDS(paste0(out.dir,'/Step-Prefilter_result_till84909.RDS'))
prefilter_df_doi <- prefilter_df %>% filter(D2 %in% disease_oi) %>% mutate(prefilter.p.adjust = p.adjust(prefilter.p, method = "bonferroni"))
pairs_to_process <- prefilter_df %>% filter(prefilter.p.adjust < 0.05) 


pairs_to_process_doi <- prefilter_df_doi %>% filter(prefilter.p.adjust < 0.05)  %>% filter((D2 %in% disease_oi))
dim(pairs_to_process_doi)

## B. Calculate RR and direction ####
SIGNIFICANCE_THRESHOLD <- 0.05 / nrow(pairs_to_process) # Bonferroni corrected threshold  
results <- list()
for(i in 1:row(pairs_to_process_doi)) { #nrow(danish_study_diseases_pairs_filtered)
  d1_code <- pairs_to_process_doi$D1[i]
  d2_code <- pairs_to_process_doi$D2[i]
  D1_name = disease_n_patient %>% filter(ICD10 == d1_code) %>% pull(ICD10_Disease)
  D2_name = disease_n_patient %>% filter(ICD10 == d2_code) %>% pull(ICD10_Disease)
  
  message('Processing disease pair: ', D1_name, ' -> ', D2_name)
  
  ### Load D1 and D2 specific data  ####
  file_d1 <- list.files(input.dir,pattern = paste0("ukb_", d1_code,"_"), full.names = TRUE)
  file_d2 <- list.files(input.dir,pattern = paste0("ukb_", d2_code,"_"), full.names = TRUE)
  
  df_d1 <- readRDS(file_d1)  %>% mutate(eid = as.character(eid))
  df_d2 <- readRDS(file_d2)  %>% mutate(eid = as.character(eid))
  
  
  ### Run RR (exposed is Patients with D1) ####
  exposed_ids <- df_d1$eid
  n_exposed <- length(exposed_ids)
  if (n_exposed < MIN_EXPOSED) {
    message("Too few exposed cases, skipping")
    next
  }
  
  exposed_data <- All_clinic %>%
    filter(eid %in% exposed_ids) %>%
    left_join(df_d1 %>% select(eid, diagnosis_time), by = "eid") %>%
    mutate(exposed = 1)
  
  # 3. Define Potential Control Pool (Individuals without D1)
  control_pool <- All_clinic %>%
    filter(!(eid %in% exposed_ids)) %>%
    mutate(exposed = 0, diagnosis_time = NA)
  
  # Combine for matching
  match_data <- bind_rows(exposed_data, control_pool)
  
  # 4. Calculate RR (Hybrid Method to handle Prevalence)
  if (n_exposed < 5000) {
    # Small Group: Use MatchIt
    match_data <- bind_rows(exposed_data, control_pool %>% mutate(exposed = 0))
    mod_match <- matchit(
      exposed ~ sex_f31_0_0 + age_at_recruitment_f21022_0_0 + 
        assessment_year + assessment_month + 
        uk_biobank_assessment_centre_f54_0_0, 
      data = match_data, 
      method = "nearest", 
      ratio = 100,
      caliper = c(age_at_recruitment_f21022_0_0 = 3, 
                  assessment_year = 3, 
                  assessment_month = 3),
      std.caliper = FALSE 
    )
    matched_controls <- match.data(mod_match) %>% filter(exposed == 0)
    
    # Count D2 in controls
    control_ids <- matched_data %>% filter(exposed == 0) %>% pull(eid)
    c_controls_total <- sum(control_ids %in% df_d2$eid)
    
    # RR = C_exposed / (Average C_controls) 
    # Since we have 100 matches per 1 exposed, the average is Total_D2_in_controls / 100
    avg_c_control <- c_controls_total / 100
  } else {
    # Large Group: Stratified Probability (Bernoulli Trial logic)
    strata_risk <- control_pool %>%
      mutate(has_d2 = ifelse(eid %in% df_d2$eid, 1, 0)) %>% select(sex_f31_0_0, max_age,has_d2) %>%
      mutate(sex_f31_0_0 = as.character(sex_f31_0_0)) %>%
      dplyr::group_by(sex_f31_0_0,max_age) %>%
      dplyr::summarise(count = n(),p_d2 = mean(has_d2)) %>%
      dplyr::ungroup()
    
    avg_c_control <- exposed_data %>%
      dplyr::group_by(sex_f31_0_0, max_age) %>%
      tally() %>%
      dplyr::left_join(strata_risk, by = c("sex_f31_0_0", "max_age")) %>%
      dplyr::mutate(exp = n * p_d2) %>%  ungroup() %>% 
      dplyr::summarise(total_exp = sum(exp, na.rm = TRUE)) %>% pull(total_exp)
  }
  
  c_exposed <- sum(df_d1$eid %in% df_d2$eid)
  rr_val <- c_exposed / avg_c_control
  
  ### Test significance ####
  # Convert avg_c_control back to a probability: p = expected_count / n_exposed
  p_expected_matched <- avg_c_control / n_exposed
  
  # Binomial test models the single comparison of the pair
  p.val <- binom.test(c_exposed_lifetime, n_exposed, p = p_expected_matched, alternative = "greater")$p.value
  
  if (p.val > SIGNIFICANCE_THRESHOLD) {
    results[[i]] <- data.frame(
      D1 = d1_code,
      D2 = d2_code,
      p_val = p_val_significance,
      RR = rr_val,
      p.value.direction = p_val_dir,
      D1_counts = n_exposed,
      D2_counts = nrow(df_d2),
      D1_name = D1_name,
      D2_name = D2_name,
      Note = "Pair not significant")
    message("Pair not significant. Skipping.")
    next
  }
  
  ### Directionality Test (Binomial) ####
  # Find cases where D1 occurred before D2 vs D2 before D1
  # This requires joining the diagnosis_time from both files
  timing_df <- inner_join(df_d1 %>% select(eid,diagnosis_time), df_d2 %>% select(eid,diagnosis_time), by = "eid", suffix = c("_d1", "_d2"))
  d1_before_d2 <- sum(timing_df$diagnosis_time_d1 < timing_df$diagnosis_time_d2)
  d2_before_d1 <- sum(timing_df$diagnosis_time_d2 < timing_df$diagnosis_time_d1)
  
  p_val_dir <- binom.test(d1_before_d2, (d1_before_d2 + d2_before_d1), p = 0.5)$p.value
  
  # Store Results
  results[[i]] <- data.frame(
    D1 = d1_code,
    D2 = d2_code,
    p_val = p_val_significance,
    RR = rr_val,
    p.value.direction = p_val_dir,
    D1_counts = n_exposed,
    D2_counts = nrow(df_d2),
    D1_name = D1_name,
    D2_name = D2_name,
    Note = NA
  )
}
saveRDS(results, file = paste0(out.dir,'/Step-RR_direction_results_doi.RDS'))
results <- readRDS(file = paste0(out.dir,'/Step-RR_direction_results_doi.RDS'))

final_trajectories <- bind_rows(results)
# Adjust p-values
final_trajectories <- final_trajectories %>%
  mutate(p.value.direction.adj = p.adjust(p.value.direction, method = "BH")) %>% 
  mutate(direction_yes_no = ifelse(p.value.direction.adj < 0.05, '1', '-1'))
final_trajectories
write_xlsx(final_trajectories, path = paste0(out.dir,'/Final_trajectories.xlsx'))


