library(dplyr)
library(stringr)
library(plyr)
library(tibble)
library(readxl)
library(dplyr)
library(purrr)
library(future.apply)

# output folder
output_folder = 'output/ICD10_level3_diseases'
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# input folder
location = 'tetralith' #omika or tetralith
input_folder = 'input'
if (location == 'mac') {
  commonfiles = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB'
}

if (location == 'tetralith') {
  setwd("/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3")
  commonfiles = '/proj/spatial_pre_mal/users/x_yelzh/Projects/Common files/UKBB' # tetralith
}

if (location == 'omika') {
  setwd( '/home/yelzh67/Projects/DiseaseTrajectory/UKBB ICD10level3')
  commonfiles = '/home/yelzh67/Projects/DiseaseTrajectory/UKBB/input'
}

##################
# #Get ICD code ####
# ICD10_coding19 <- read.csv2(paste0(input_folder,'/UKBB ICD10 code dic/coding19.tsv'),sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# head(ICD10_coding19)
# 
# # Select ICD10 code level3 (select code with 3 characters in coding column)
# ICD10_coding19_level3 <- ICD10_coding19 %>%
#   filter(nchar(coding) == 3) %>%
#   select(coding, meaning) %>%
#   distinct()
# dim(ICD10_coding19_level3) # total of 2066 diseases at level 3 which will have 4,266,290 pairs (consider direction)
# 
# #check the names of disease_list, replace to '_' if there are ' ', ',' , '/' , '(',')',  in the names
# ICD10_coding19_level3$meaning <- ICD10_coding19_level3$meaning |>
#   gsub("[ ,/&]", "_", x = _) |>     # replace space, comma, slash, ampersand with "_"
#   gsub("[()'\\.;]", "", x = _) |>   # remove parentheses, apostrophe, dot, semicolon
#   gsub("_+", "_", x = _)            # collapse multiple underscores into one
# write.csv2(ICD10_coding19_level3, file = paste0(output_folder,"/icd10_level3_disease_names.csv"), row.names = FALSE)
ICD10_coding19_level3 = read.csv2(paste0(output_folder,"/icd10_level3_disease_names.csv"), header = TRUE, stringsAsFactors = FALSE) 
# # create a named list of diseases with their corresponding ICD10 codes
# disease_list <- split(ICD10_coding19_level3$coding,ICD10_coding19_level3$meaning)
# 
# save(disease_list, file = paste0(output_folder,"/icd10_level3_disease_list.RData"))

#Load the disease list
load(paste0(output_folder,"/icd10_level3_disease_list.RData"))

# transfer file to server####
# scp -r "/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB"  x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/Common files'
# scp -r "/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3"  x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/'
# scp -r "/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases"  ssh yelzh67@10.237.74.254:'/home/yelzh67/Projects/DiseaseTrajectory'
# transfer output from tetralith to mac to omika
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output'
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'
# scp -r "/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases"  ssh yelzh67@10.237.74.254:'/home/yelzh67/Projects/DiseaseTrajectory/UKBB ICD10level3/output/'
# transfer output from omika to mac to tetralith
# scp -r yelzh67@10.237.74.254:'/home/yelzh67/Projects/DiseaseTrajectory/UKBB ICD10level3/output/ICD10_level3_diseases'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output'
# scp -r "/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases"  x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output'
# transfer output between omika and tetralith
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases' yelzh67@10.237.74.254:'/home/yelzh67/Projects/DiseaseTrajectory/UKBB ICD10level3/output'
# scp -r yelzh67@10.237.74.254:'/home/yelzh67/Projects/DiseaseTrajectory/UKBB ICD10level3/output/ICD10_level3_diseases'  x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output'

################## Read data ###################################
## Choose either genomics, proteomics and metabolomics
omics_path = paste(commonfiles, '/Olink_proteomics_data_2ndPhase_transposed_decoded2UNIportID.txt', sep = '')
output_file = 'prot' #is used in saving results

#load all UKBB patients
load(paste(commonfiles, '/ukb672643.rda', sep = ''))
dim(my_ukb_data) #[1] 502369  18518

# if(output_file == 'prot'){
#   PID_omics = read.table(file = omics_path, sep='\t', header = TRUE, fill = TRUE, row.names = 'PID')
#   dim(PID_omics) #[1] 53073  2923
#   PID_omics = rownames(PID_omics)
# }

#clinical with columns only relevant for finding healthy controls
# diagnosis_columns <- paste(commonfiles, '/Decoded_ukb_Matrix_Firoj_ControlHealthyICD9and10codes_only.csv', sep = '')
# ukb672643_diagnosis_columns <- read.csv(diagnosis_columns, header = TRUE, stringsAsFactors = FALSE)
# save(ukb672643_diagnosis_columns, file = paste0(commonfiles, '/ukb672643_diagnosis_columns.rda'))
load(paste0(commonfiles, '/ukb672643_diagnosis_columns.rda'))

# A list of patients that withdrawn from the UKBB
withdrawn_patients = read.table(paste(commonfiles, '/UKBBwithdraw.csv', sep = ''))
my_ukb_data = my_ukb_data[!my_ukb_data$eid %in% withdrawn_patients$V1,] 
dim(my_ukb_data) #  [1] 502359  18518
ukb672643_diagnosis_columns = ukb672643_diagnosis_columns[!ukb672643_diagnosis_columns$eid %in% withdrawn_patients$V1,]
dim(ukb672643_diagnosis_columns) #  [1] 502359  736
# PID_omics = PID_omics[!PID_omics %in% withdrawn_patients$V1]
# length(PID_omics) #

################## Define varibales of interests ##############################
clinical_variables = c('eid','year_of_birth_f34_0_0', 'age_at_recruitment_f21022_0_0','date_of_attending_assessment_centre_f53_0_0',
                       'sex_f31_0_0', 'ethnic_background_f21000_0_0','uk_biobank_assessment_centre_f54_0_0',
                       'hand_grip_strength_left_f46_0_0', 'hand_grip_strength_right_f47_0_0', 
                       'smoking_status_f20116_0_0', 'alcohol_drinker_status_f20117_0_0',"alcohol_intake_frequency_f1558_0_0"  ,              
                       'diabetes_diagnosed_by_doctor_f2443_0_0', 'diastolic_blood_pressure_automated_reading_f4079_0_0',
                       'systolic_blood_pressure_automated_reading_f4080_0_0', 'body_mass_index_bmi_f21001_0_0', 'glucose_f30740_0_0',
                       'haemoglobin_concentration_f30020_0_0',
                       'haematocrit_percentage_f30030_0_0', 'mean_corpuscular_volume_f30040_0_0', 'mean_corpuscular_haemoglobin_f30050_0_0',
                       'red_blood_cell_erythrocyte_distribution_width_f30070_0_0', 'platelet_crit_f30090_0_0', 'mean_platelet_thrombocyte_volume_f30100_0_0',
                       'platelet_distribution_width_f30110_0_0', 'platelet_distribution_width_f30110_0_0', 'lymphocyte_percentage_f30180_0_0', 
                       'monocyte_percentage_f30190_0_0', 'neutrophill_percentage_f30200_0_0', 'basophill_percentage_f30220_0_0', 'reticulocyte_percentage_f30240_0_0',
                       'mean_reticulocyte_volume_f30260_0_0', 'immature_reticulocyte_fraction_f30280_0_0', 'high_light_scatter_reticulocyte_percentage_f30290_0_0',
                       'creatinine_enzymatic_in_urine_f30510_0_0', 'potassium_in_urine_f30520_0_0', 'sodium_in_urine_f30530_0_0', 'alkaline_phosphatase_f30610_0_0', 
                       'alanine_aminotransferase_f30620_0_0', 'apolipoprotein_b_f30640_0_0', 'aspartate_aminotransferase_f30650_0_0', 'urea_f30670_0_0',
                       'cholesterol_f30690_0_0', 'creatinine_f30700_0_0', 'creactive_protein_f30710_0_0', 
                       'gamma_glutamyltransferase_f30730_0_0', 'glycated_haemoglobin_hba1c_f30750_0_0', 'igf1_f30770_0_0', 'ldl_direct_f30780_0_0',
                       'total_bilirubin_f30840_0_0',  'triglycerides_f30870_0_0', 'urate_f30880_0_0',
                       "date_of_death_f40000_0_0","date_of_death_f40000_1_0","underlying_primary_cause_of_death_icd10_f40001_0_0","age_at_death_f40007_0_0","age_at_death_f40007_1_0",
                       "reason_lost_to_followup_f190_0_0","date_lost_to_followup_f191_0_0")

save(clinical_variables, file = paste0(output_folder,'/ukb_clinic.rda'))
# my_ukb_data[my_ukb_data$eid == '1000829',grep("^diagnoses.*icd10", colnames(my_ukb_data))]


################## Identify patients ###################################
# extract patients with certain ICD code ####
ukb_icd_subset_by_ICD = function (data, icd.code, icd.version = 10) 
{
  ukb_case <- data %>% dplyr::select(matches(paste("^diagnoses.*icd", 
                                                   icd.version, sep = ""))) %>% 
    purrr::map_df(~grepl(icd.code, ., perl = TRUE)) %>% rowSums() > 0
  data_subset <- data[ukb_case,]
  return(data_subset)
}

# variables that relates to diagnosis
icd.version = '10'

# Repeat for all diseases
for(j in 1:length(disease_list)){ #length(disease_list)
  disease = names(disease_list[j] )
  icd.code.list = disease_list[[j]]
  message(paste0("Processing disease No ",j,": ", disease, " with ICD codes: ", paste(icd.code.list, collapse = ", ")))
  
  if(file.exists(paste0(output_folder, '/ukb_', disease, '_subset_clinic.rda'))){
    message(paste0("File for ", disease, " already exists. Skipping..."))
    next
  }
  
  ukb_disease_subset_list = list()
  ukb_case_list = list()
  for(icd.code in icd.code.list){
    # all patients with a given icd.code
    patients_with_disease = ukb_icd_subset_by_ICD(my_ukb_data, icd.version = icd.version, icd.code = icd.code)
    
    # the below is to remember which columns corresponded to the disease and is later used for creating ToBeSick group
    ukb_case_subset <- patients_with_disease %>% dplyr::select(matches("40006|41270")) %>%
      purrr::map_df(~grepl(icd.code, ., perl = TRUE))
    
    if(dim(ukb_case_subset)[1] == 0){
      print(paste0("no patients with disease ICD code ", icd.code))
      next
    }
    ukb_case_subset = cbind(eid = patients_with_disease$eid,ukb_case_subset)
    
    ukb_case_list = append(ukb_case_list, list(ukb_case_subset))
    
  }

  ukb_case = bind_rows(ukb_case_list)
  
  if(dim(ukb_case)[1] == 0){
    print("no patients with disease")
    next
  }
  
  # There might be duplicates in ukb_case (if for example patient has both M05 and M06 in RA)
  # The below code takes only unique values and for each diagnosis returns TRUE if it belongs to either M05 or M06
  coln = colnames(ukb_case)
  ukb_case = ddply(ukb_case, .(eid), function(x) return(as.logical(colSums(x[,2:dim(ukb_case)[2]]))))
  #ukb_case <- ukb_case %>%  group_by(eid) %>%  summarise(across(everything(), ~ any(.)), .groups = "drop")
  colnames(ukb_case) = coln
  rownames(ukb_case) = ukb_case$eid
  ukb_case = ukb_case[,-1]
  
  ukb_disease_subset = my_ukb_data[my_ukb_data$eid %in% rownames(ukb_case),]
  
  ## Identify DISEASE diagnosis date  ###################################
  # for each patient, we find all disease related columns, then we map them into corresponding column names for diagnosis time,
  # then remove all columns with secondary disease and select the earliest diagnosis time
  earliest_diagnosis = data.frame()
  for(i in 1:dim(ukb_disease_subset)[1]){
    diagnoses = colnames(ukb_case)[unlist(ukb_case[i,])]
    # if diagnoses contains '41270' then run below code diagnoses = gsub("type_of_cancer_icd10_f40006", "date_of_cancer_diagnosis_f40005", diagnoses) # This is because of the codes for diagnosis and time of diagnosis is different
    if (any(grepl("41270", diagnoses))) {
      diagnoses <- gsub("diagnoses_icd10_f41270", "date_of_first_inpatient_diagnosis_icd10_f41280", diagnoses)
    } else if (any(grepl("41271", diagnoses)))  {
      diagnoses <- gsub("diagnoses_icd9_f41271", "date_of_first_inpatient_diagnosis_icd9_f41281", diagnoses)
    }
    if (any(grepl("40006", diagnoses))) {
      diagnoses <- gsub("type_of_cancer_icd10_f40006", "date_of_cancer_diagnosis_f40005", diagnoses)
    } 
    
    diagnoses_times = ukb_disease_subset[i,diagnoses]
    
    if(length(diagnoses_times) == 1){
      earliest_diagnosis[i,'eid'] = ukb_disease_subset[i,'eid']
      earliest_diagnosis[i,'disease'] = disease
      earliest_diagnosis[i,'diagnosis_time'] = c(diagnoses_times)
    }else{
      earliest_diagnosis[i,'eid'] = ukb_disease_subset[i,'eid']
      earliest_diagnosis[i,'disease'] = disease
      earliest_diagnosis[i,'diagnosis_time'] = min(as.Date(as.vector(t(diagnoses_times)), format = "%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  ## Extract clinical columns
  ukb_disease_subset = merge(earliest_diagnosis,my_ukb_data[,clinical_variables], by = 'eid', all.x = TRUE)
  ukb_disease_subset = ukb_disease_subset %>% mutate(prev_inci = as.numeric(difftime(diagnosis_time, date_of_attending_assessment_centre_f53_0_0, units = "days"))<=0, 'prev','inci') 
  print(paste0(disease, ' ICD10: ',paste(icd.code.list, collapse = ", "),' Number of patients: ', dim(ukb_disease_subset)[1]))
  # #write.csv2(ukb_disease_subset, file = paste0(output_folder, '/ukb_', disease, '_subset_clinic.csv'), row.names = FALSE)
  # saveRDS(ukb_disease_subset[ukb_disease_subset$prev_inci == 'prev',], file = paste0(output_folder, '/ukb_', disease, '_prev_subset_clinic.rda'))
  # saveRDS(ukb_disease_subset[ukb_disease_subset$prev_inci == 'inci',], file = paste0(output_folder, '/ukb_', disease, '_inci_subset_clinic.rda'))
  saveRDS(ukb_disease_subset, file = paste0(output_folder, '/ukb_', disease, '_subset_clinic.rda'))
}

################## Identify HC that without any diagnosis ###################################
baseline = 'HC'
if(baseline == 'HC'){
  rows_with_all_na = ukb672643_diagnosis_columns[rowSums(is.na(ukb672643_diagnosis_columns[, -1])) == (ncol(ukb672643_diagnosis_columns) - 1), ]
  ukb_healthy_controls = my_ukb_data[my_ukb_data$eid %in% rows_with_all_na$eid,clinical_variables]
  saveRDS(ukb_healthy_controls, file = paste(output_folder, '/ukb_healthy_controls.rda', sep = ''))
}
 

