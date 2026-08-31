# In this version, I used the same method for Danish data mentioned in https://doi.org/10.1038/s41467-020-18682-4
# But instead of using their results, I run the whole pipe line again to get RR and directionality p-values in UKBB for trajectories towards to common cancers
# This is a updated version on 20260620 that cleaned for below document https://kise.sharepoint.com/:w:/r/teams/GRP_MDT/Delade%20dokument/EARLY%20T/Prong%201%20Texts.docx?d=w31e980be674e4ba2b90c00eec4e08f62&csf=1&web=1&e=nUmt0M
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
out.dir1 = "output/ICD10_level3_diseases/trajectory_vali"
out.dir = "output/ICD10_level3_diseases/trajectory_vali_20260623"
dir.create(out.dir, recursive = T,showWarnings = FALSE)
input.dir = 'output/ICD10_level3_diseases'
comm_file.dir = '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/Common files/UKBB'

location = 'tetralith'
if (location == 'tetralith') {
  setwd("/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3")
  comm_file.dir = '/proj/spatial_pre_mal/users/x_yelzh/Projects/Common files/UKBB' # tetralith
}

################## 
#1 Load diseases in UKBB ####
disease_n_patient <- readRDS(file = paste0(input.dir,'/disease_n_patient.RDS'))
disease_n_patient_Prot <- readRDS(file = paste0(input.dir,'/disease_n_patient_proteomics.RDS'))
sex_specific_diseases <- readRDS(file = paste0(input.dir,'/sex_specific_diseases.rds'))

#2 Load UKBB clinic data ####
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

#3 Load pre-filter results  ####
prefilter_df <- readRDS(file = paste0(out.dir1,'/Step-Prefilter_result.RDS'))

#4 Subset trajectories to end diseases of interest ####
disease_oi <- c("C50",# Breast C
                "C61",# Prostate C
                "C18",# Colon C
                "C19",#Malignant_neoplasm_of_rectosigmoid_junction
                "C20",#Malignant_neoplasm_of_rectum"
                "C43",#C43_Malignant_melanoma_of_skin
                "C53",#C53_Malignant_neoplasm_of_cervix_uteri
                "C54",#Malignant_neoplasm_of_corpus_uteri
                "C55",#Malignant_neoplasm_of_uterus_part_unspecified
                "C62",#Malignant_neoplasm_of_testis
                "C56",#Malignant_neoplasm_of_ovary
                "C67",#Malignant_neoplasm_of_bladder
                "C64",#Malignant_neoplasm_of_kidney_except_renal_pelvis
                "C34",#Malignant_neoplasm_of_bronchus_and_lung
                "C25", #Malignant_neoplasm_of_pancreas
                "C22" #C22 Malignant neoplasm of liver and intrahepatic bile ducts
) 

prefilter_df_doi <- prefilter_df %>% filter(D2 %in% disease_oi) %>% 
  mutate(prefilter.p.adjust = p.adjust(prefilter.p, method = "fdr")) #19184
# save it to excel file
write.csv(prefilter_df_doi,paste0(out.dir,'/prefilter_df_doi_non_filted.csv'))


pairs_to_process <- prefilter_df_doi %>%  
  mutate(prefilter.p.adjust = p.adjust(prefilter.p, method = "fdr")) %>% 
  filter(prefilter.p.adjust < 0.05)  #247
dim(pairs_to_process)

#5 Calculate RR and direction ####
results <- list()
for(i in 1:nrow(pairs_to_process)) { 
  d1_code <- pairs_to_process$D1[i]
  d2_code <- pairs_to_process$D2[i]
  
  # Get names 
  D1_name <- disease_n_patient %>% filter(ICD10 == d1_code) %>% pull(ICD10_Disease)
  D2_name <- disease_n_patient %>% filter(ICD10 == d2_code) %>% pull(ICD10_Disease)
  
  message('Processing disease pair ',i, ':', D1_name, ' -> ', D2_name)
  
  ### 1. Load Data ####
  file_d1 <- list.files(input.dir, pattern = paste0("ukb_", d1_code,"_"), full.names = TRUE)
  file_d2 <- list.files(input.dir, pattern = paste0("ukb_", d2_code,"_"), full.names = TRUE)
  
  df_d1 <- readRDS(file_d1[1]) %>% mutate(eid = as.character(eid))
  df_d2 <- readRDS(file_d2[1]) %>% mutate(eid = as.character(eid))
  
  # Define exposed (Patients with D1)
  exposed_data <- All_clinic %>% 
    filter(eid %in% df_d1$eid) %>% left_join(
      df_d1 %>% select(eid, diagnosis_time), by = "eid"
    ) %>%
    mutate(exposed = 1)
  
  n_exposed <- nrow(exposed_data)
  
  ### 2. Matching (N=10 as requested) ####
  # Define control pool (Individuals without D1)
  control_pool <- All_clinic %>% 
    filter(!(eid %in% df_d1$eid)) %>%
    mutate( diagnosis_time = NA, exposed = 0)
  
  match_data_combined <- bind_rows(exposed_data, control_pool)
  
  # Matching logic using your variables
  mod_match <- tryCatch({
    matchit(
      exposed ~ sex_f31_0_0 + max_age + assessment_year + 
        uk_biobank_assessment_centre_f54_0_0, 
      data = match_data_combined, 
      method = "nearest", 
      ratio = 10,  # N=10 matching
      caliper = c(assessment_year = 3),
      std.caliper = FALSE 
    )
  }, error = function(e) return(NULL))
  
  if (is.null(mod_match)) {
    message("Matching failed for this pair. Skipping.")
    next
  }
  
  # Extract matched controls
  matched_data <- match.data(mod_match)
  matched_controls <- matched_data %>% filter(exposed == 0)
  ratio_matched <- n_exposed/nrow(matched_controls) 
  
  ### 3. Calculate RR ####
  # Join to find any sequence where D1 occurred before D2
  timing_df <- inner_join(
    df_d1 %>% select(eid, diagnosis_time_d1 = diagnosis_time), 
    df_d2 %>% select(eid, diagnosis_time_d2 = diagnosis_time), 
    by = "eid"
  )
  
  # C_exposed: Cases of D2 after D1 
  c_exposed <- timing_df %>% 
    filter(diagnosis_time_d2 > diagnosis_time_d1) %>% 
    nrow()
  
  # C_controls_total: Total cases of D2 found in the matched control group
  c_controls_total <- sum(matched_controls$eid %in% df_d2$eid)
  
  # avg_c_control 
  avg_c_control <- c_controls_total * ratio_matched
  
  # Calculate RR (avoid division by zero if controls have 0 cases)
  rr_val <- ifelse(avg_c_control > 0, c_exposed / avg_c_control, NA)
  
  ### 4. Test Significance (Binomial Distribution) ####
  # Expected probability p = avg_expected_count / n_exposed
  p_expected <- avg_c_control / n_exposed
  
  # Paper logic: binomial test models the comparison
  p_val_significance <- binom.test(c_exposed, n_exposed, p = p_expected, alternative = "greater")$p.value
  
  ### 5. Directionality Test ####
  # Find temporal sequence for the directionality step
  d1_before_d2 <- sum(timing_df$diagnosis_time_d1 < timing_df$diagnosis_time_d2, na.rm = TRUE)
  d2_before_d1 <- sum(timing_df$diagnosis_time_d2 < timing_df$diagnosis_time_d1, na.rm = TRUE)
  
  # Paper: binomial distribution for directionality
  p_val_dir <- binom.test(d1_before_d2, (d1_before_d2 + d2_before_d1), p = 0.5, alternative = "greater")$p.value
  
  # Store Results
  results[[i]] <- data.frame(
    D1 = d1_code,
    D2 = d2_code,
    D1_name = D1_name,
    D2_name = D2_name,
    RR = rr_val,
    p_val_RR = p_val_significance,
    p_val_direction = p_val_dir,
    D1_counts = n_exposed,
    D2_counts = nrow(df_d2),
    D1_before_D2_counts = d1_before_d2,
    d2_before_D1_counts = d2_before_d1,
    Delta_D2_D1_diagnosis_time = mean(timing_df$diagnosis_time_d2 - timing_df$diagnosis_time_d1, na.rm = TRUE)
    # Direction_Significant = ifelse(p_val_dir < 0.05, "Yes", "No")
  )
}
saveRDS(results, file = paste0(out.dir,'/Step-RR_direction_results_doi_list.RDS'))
# Bind into a single dataframe
final_results_df <- bind_rows(results)
saveRDS(final_results_df, file = paste0(out.dir,'/Step-RR_direction_results_doi.RDS'))

final_results_df <- readRDS(file = paste0(out.dir,'/Step-RR_direction_results_doi.RDS'))
write_xlsx(final_results_df, path = paste0(out.dir,'/Final_trajectories_UKBB_doi_not-filtered.xlsx'))
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'
# scp -r '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali' x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'  
# scp -r x_yelzh@tetralith.nsc.liu.se:'/proj/spatial_pre_mal/users/x_yelzh/Projects/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases/trajectory_vali_20260623'  '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/UKBB ICD10level3/output/ICD10_level3_diseases'

final_trajectories <- final_results_df %>% mutate(
  p.value.RR.adj = p.adjust(p_val_RR, method = "fdr"),
  RR_significant = ifelse(p.value.RR.adj < 0.05, '1', '-1')) %>%
  filter(RR_significant == '1') 

dim(final_trajectories)

final_trajectories <- final_results_df %>% mutate(
  p.value.RR.adj = p.adjust(p_val_RR, method = "fdr"),
  RR_significant = ifelse(p.value.RR.adj < 0.05, '1', '-1')) %>%
  filter(RR_significant == '1') %>% 
  mutate(p_val_direction.adj = p.adjust(p_val_direction, method = "fdr"),
  direction_yes_no = ifelse(p_val_direction.adj < 0.05, '1', '-1')) %>% 
  filter(direction_yes_no == '1') %>%
  arrange(D2,desc(RR))
write_xlsx(final_trajectories, path = paste0(out.dir,'/Final_trajectories_UKBB_doi.xlsx'))

final_trajectories %>% dim

# Make sanky plot from D1 to D2####
library(networkD3)
sankey_data <- final_trajectories %>%
  select(D1_name, D2_name, RR) %>%
  mutate(RR_scaled = scales::rescale(RR, to = c(1, 10))) # Scale RR for better visualization
# Create nodes
nodes <- data.frame(name = unique(c(sankey_data$D1_name, sankey_data$D2_name)))
# Create links 
links <- sankey_data %>%
  mutate(
    source = match(D1_name, nodes$name) - 1,
    target = match(D2_name, nodes$name) - 1,
    # 'value' controls line thickness. Use your scaled RR here.
    value = RR_scaled, 
    # Add the actual RR for the tooltip display
    RR_display = round(RR, 2) 
  ) %>%
  select(source, target, value, RR_display)

# 3. Plot Sankey diagram
# We set 'Value = "RR_display"' if you want the tooltip to show the RR.
# Note: In Sankey plots, the line width is tied to this value. 
p <- sankeyNetwork(
  Links = links, 
  Nodes = nodes,
  Source = "source", 
  Target = "target",
  Value = "RR_display",  # This ensures the RR value is what users see
  NodeID = "name",
  units = "Relative Risk", # This text appears after the number in the tooltip
  fontSize = 8, 
  nodeWidth = 50,
  nodePadding = 6,         # ADDED: Adds vertical space between nodes
  height = 900,            # ADDED: Increases the overall height of the plot
  width = 1000,
  sinksRight = FALSE      # Often looks better for trajectory flow
)

# 4. Save the plot
htmlwidgets::saveWidget(p, file = paste0(out.dir, '/Sankey_plot_doi.html'))

# Bubble plot ####
library(ggplot2)
library(ggrepel)
library(igraph)
library(cowplot)
final_trajectories_plot_df = final_trajectories %>%
  mutate(
    type = ifelse(grepl("^C", D1), "Cancer-to-Cancer", 
                  ifelse(grepl("^D", D1), "Benign-to-Cancer", "Non-Neoplastic-Disease-to-Cancer")),
    label = paste0(D1, "→", D2)
  )


p = final_trajectories_plot_df  %>%
  ggplot(aes(x = Delta_D2_D1_diagnosis_time/365,  
             y = log2(RR), size = log2(RR), color = type, label = label)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_text_repel(size = 4, max.overlaps = 10) +
  scale_size_continuous(range = c(3, 14), guide = "none") +
  scale_color_manual(values = c("Benign-to-Cancer" = "#378ADD",
                                "Cancer-to-Cancer"    = "#1D9E75", 
                                "Non-Neoplastic-Disease-to-Cancer"    = "orange")) +
  labs(x = "Mean time D2 - D1 (years)", y = "Relative Risk (log2)",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank())

# 2. Extract the legend from the base plot
my_legend <- get_legend(p)

# 3. Strip the legend from the base plot so ggMarginal doesn't overlay it
p_no_legend <- p + theme(legend.position = "none")

# 4. Generate the marginal plots on the clean, legend-less plot
p_marginal = ggMarginal(
  p_no_legend, 
  type = "density", 
  groupColour = TRUE, 
  groupFill = TRUE,   
  alpha = 0.3         
)

final_plot <- plot_grid(p_marginal, my_legend, ncol = 2, rel_widths = c(1, 0.25))

# 6. Save the plot
pdf(file = paste0(out.dir,"/UKBB_doi_sig_RR_and_time_diff_Bubble_Type.pdf"), width = 12, height = 8)
print(final_plot)
dev.off()

# Bubble plot - color by DIS2 ####
library(cowplot)
# 1. Create your base plot normally 
p = final_trajectories_plot_df %>% filter(type == "Disease1 → Cancer") %>%
  ggplot(aes(x = Delta_D2_D1_diagnosis_time/365,  
             y = RR, size = RR, color = D2_name, label = label)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_text_repel(size = 4, max.overlaps = 10) +
  scale_y_log10() +
  scale_size_continuous(range = c(3, 14), guide = "none") +
  labs(x = "Mean time D2 - D1 (years)", y = "Relative Risk",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank())

my_legend <- get_legend(p)
p_no_legend <- p + theme(legend.position = "none")
p_marginal = ggMarginal(
  p_no_legend, 
  type = "density", 
  groupColour = TRUE, 
  groupFill = TRUE,   
  alpha = 0.3         
)
 
final_plot <- plot_grid(p_marginal, my_legend, ncol = 2, rel_widths = c(1, 0.4))

# 6. Save the plot 
# NOTE: I increased the width from 12 to 14 to give the long legend text breathing room
pdf(file = paste0(out.dir,"/UKBB_doi_sig_RR_and_time_diff_Bubble_D2.pdf"), width = 14, height = 8)
print(final_plot)
dev.off()