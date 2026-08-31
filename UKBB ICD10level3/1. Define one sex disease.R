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

# setup folder ####
out.dir = 'output/ICD10_level3_diseases'
dir.create(out.dir, recursive = T,showWarnings = FALSE)
input.dir = 'output/ICD10_level3_diseases'
comm_file.dir = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB'
disease_files <- list.files(out.dir, pattern = "_subset_clinic.rda", full.names = TRUE)

# 1. Cleaning - misclassified diseases ####
library(dplyr)
library(stringr)

# Define the subfolder for "original" (mismatched) files
archive_subfolder <- "original_with_mismatches"
dir.create(file.path(out.dir, archive_subfolder), showWarnings = FALSE)

clean_sex_mismatch <- function(file_path, expected_sex) {
  df <- readRDS(file_path)
  df_name <- basename(file_path)
  
  # Identify misclassified patients 
  mismatched <- df %>% filter(sex_f31_0_0 != expected_sex)
  
  if (nrow(mismatched) > 0) {
    message(paste("Mismatches found in:", df_name, "| Count:", nrow(mismatched)))
    
    # 1. Document the mismatches for your log
    temp_log <- data.frame(
      file = df_name,
      expected_sex = expected_sex,
      mismatched_count = nrow(mismatched),
      mismatched_ids = paste(mismatched$eid, collapse = ";")
    )
    
    # 2. Define the path for the ARCHIVE (moving the old file)
    archive_path <- file.path(dirname(file_path), archive_subfolder, df_name)
    
    # 3. Move the old file to the subfolder
    file.rename(from = file_path, to = archive_path)
    
    # 4. Filter the data
    df_cleaned <- df %>% filter(sex_f31_0_0 == expected_sex)
    
    # 5. Save the NEW (cleaned) file to the ORIGINAL location
    saveRDS(df_cleaned, file_path)
    
    return(temp_log)
  }
  
  return(NULL)
}


# define keywords
male_keywords <- c("prostate", "testis", "testicular", "penis", "seminal", 
                   "scrotum", "epididymis", "prepuce", "azoospermia", "erectile","_male_genital","_male")

female_keywords <- c("cervix", "uteri", "uterus", "ovary", "ovarian", 
                     "vulva", "vagina", "salpingitis", "endometriosis", 
                     "menstruation", "menopause", "pregnan", "obstetric", 
                     "childbirth", "abortion", "leiomyoma", "endometrial", 
                     "puerperium", "placenta","_female_genital","_female")

mismatch_log <- data.frame()
for (file in disease_files) {
  file_lower <- tolower(basename(file))
  target_sex <- NULL
  
  # Determine expected sex based on filename
  if (any(sapply(male_keywords, grepl, x = file_lower))) {
    target_sex <- 'Male' # Male
  } else if (any(sapply(female_keywords, grepl, x = file_lower))) {
    target_sex <- 'Female'
  }
  
  # If a sex-specific keyword was found, run the cleaning function
  if (!is.null(target_sex)) {
    result <- clean_sex_mismatch(file, target_sex)
    if (!is.null(result)) {
      mismatch_log <- rbind(mismatch_log, result)
    }
  }
}

# Save the log of misclassified patients
write.csv(mismatch_log, paste0(out.dir, "/original_with_mismatches/misclassified_patients_log.csv"), row.names = FALSE)


#2. Identify sex specific diseases ####
sex_specific_diseases <- data.frame(ICD10_Disease = character(), ICD10 = character(),sex_type = character(), Number_pt = integer(), stringsAsFactors = FALSE)
# read disease files, if it only contains one sex, then save it to sex_specific_diseases
for (file in disease_files) {
  df = readRDS(file)
  df_name = gsub("_subset_clinic.rda", "", basename(file)) %>% gsub("ukb_", "", .)
  ICD10_code = str_split(df_name, "_")[[1]][1]
  table(df$sex_f31_0_0)
  if(length(unique(df$sex_f31_0_0)) >1) {next
  } else {sex_specific_diseases <- rbind(sex_specific_diseases, data.frame(ICD10_Disease = df_name, ICD10 = ICD10_code,
                                                                           sex_type = unique(df$sex_f31_0_0),Number_pt = nrow(df)))}
}
sex_specific_diseases = sex_specific_diseases %>% filter(Number_pt >= 20)
table(sex_specific_diseases$sex_type)
write.csv(sex_specific_diseases, paste0(out.dir, "/sex_specific_diseases.csv"), row.names = FALSE)
saveRDS(sex_specific_diseases,file = paste0(out.dir,"/sex_specific_diseases.rds"))

