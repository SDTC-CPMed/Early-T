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
library(circlize)
library(viridis)
library(readxl)
library(writexl)
library(stringr)
library(data.table)
library(glmnet)
library(reticulate)
library(VIM)
library(purrr)
library(broom)
library(pROC)
library(tidyverse)
library(foreach)
library(doParallel)
rm(list=ls())

# output folder ####
out.dir = 'output/ICD10_level3_diseases'
dir.create(out.dir, recursive = T,showWarnings = FALSE)
input.dir = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'
comm_file.dir = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB'

# Load disease list and proteomics data ####
icd10_level3_disease_names <- read.csv2(paste0(out.dir,'/icd10_level3_disease_names.csv'),sep = ';', header = TRUE, stringsAsFactors = FALSE)
head(icd10_level3_disease_names)

load(paste0(out.dir,"/icd10_level3_disease_list.RData"))
length(disease_list)  # 2066
UKBB_Proteomics = read.table(file = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/UKBB/input/Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt', sep='\t', header = TRUE, fill = TRUE, row.names = 'PID')

# Identify number of patients per disease ####
disease_files <- list.files(out.dir, pattern = "_subset_clinic.rda", full.names = TRUE)
disease_n_patient <- data.frame(ICD10_Disease = character(), ICD10 = character(),Number_pt = integer(), stringsAsFactors = FALSE)
disease_n_patient_Prot <- data.frame(ICD10_Disease = character(), ICD10 = character(), Number_pt = integer(), stringsAsFactors = FALSE)
for (file in disease_files) {
  df = readRDS(file)
  df_Prot = df %>% filter(eid %in% rownames(UKBB_Proteomics))
  
  df_name = gsub("_subset_clinic.rda", "", basename(file)) %>% gsub("ukb_", "", .)
  ICD10_code = disease_list[[df_name]] 
  disease_n_patient <- rbind(disease_n_patient, data.frame(ICD10_Disease = df_name, ICD10 = ICD10_code, Number_pt = nrow(df)))
  disease_n_patient_Prot <- rbind(disease_n_patient_Prot, data.frame(ICD10_Disease = df_name, ICD10 = ICD10_code, Number_pt = nrow(df_Prot)))
}


write.table(disease_n_patient, file = paste0(out.dir,'/number_of_patients.txt'), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(disease_n_patient, file = paste0(out.dir,'/disease_n_patient.RDS'))
disease_n_patient <- readRDS(file = paste0(out.dir,'/disease_n_patient.RDS'))

write.table(disease_n_patient_Prot, file = paste0(out.dir,'/number_of_patients_proteomics.txt'), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(disease_n_patient_Prot, file = paste0(out.dir,'/disease_n_patient_proteomics.RDS'))
disease_n_patient_Prot <- readRDS(file = paste0(out.dir,'/disease_n_patient_proteomics.RDS'))


disease_n_patient_filtered = disease_n_patient %>% filter(Number_pt >= 100)
disease_n_patient_Prot_filtered = disease_n_patient_Prot %>% filter(Number_pt >= 40)




