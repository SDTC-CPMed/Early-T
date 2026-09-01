
#This version will first identify all pathway files
# Remove duplicate files by keeping the newer one
# Then save all the pathway files in a new folder for further analysis.
# manual exam file names to unify them to make sure they can be run in pathway comparisons

library(reticulate)
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)
library(readxl)
library(dplyr)
library(circlize)
# library(clusterProfiler)  
library(DOSE)
# library(topGO) 
# remove all variables
rm(list = ls())

setwd('/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/Prong2')
inputdir <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Delade dokument - GRP_MDT_CPMed/EARLY T/Datasets' 
outdir <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/Prong2/Prong2A_Final_clean_version_from_YZ/output/Cleaned IPA results'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Identify files that start with 'Pathway_DEGs_toplogfc_' in the inputdir including subfolders ####
pathway_files <- list.files(inputdir, pattern = "^Pathway_DEGs_.*\\.txt$", full.names = TRUE, recursive = TRUE)

# Remove duplicate files by keeping the newer one and save to outdir####
for (file in pathway_files) {
  file_name <- basename(file)
  output_file <- file.path(outdir, file_name)
  
  if (file.exists(output_file)) {
    # If the file already exists, compare the creation times
    existing_time <- file.info(output_file)$ctime
    new_time <- file.info(file)$ctime
    
    if (new_time > existing_time) {
      # If the new file is newer, copy it to the output directory
      file.copy(file, output_file, overwrite = TRUE)
      message("Updated: ", file_name)
    } else {
      message("Skipped (existing file is newer): ", file_name)
    }
  } else {
    # If the file does not exist, copy it to the output directory
    file.copy(file, output_file)
    message("Copied: ", file_name)
  }
}

# Rename files in outdir to make them more unified ####
# for file names contains "_toplogfc_morethan", replace it with "_logFC"
final_files <- list.files(outdir, pattern = "^Pathway_DEGs_.*\\.txt$", full.names = TRUE)
for (file in final_files) {
  file_name <- basename(file)
  new_file_name <- gsub("_toplogfc_morethan", "_logFC", file_name)
  new_file_path <- file.path(outdir, new_file_name)

  if (file.exists(new_file_path)) {
    message("File already exists with the new name: ", new_file_name)
  } else {
    file.rename(file, new_file_path)
    message("Renamed: ", file_name, " to ", new_file_name)
  }
}

# Check the naming pattern####
# if split by '_' , check how many elements per file name in the outdir, highlight those that have more than 8 elements, print the filename
final_files <- list.files(outdir, pattern = "^Pathway_DEGs_.*\\.txt$", full.names = TRUE)
for (file in final_files) {
  file_name = basename(file)
  elements = strsplit(file_name, "_")[[1]]
  
  if (length(elements) > 8) {
    message("File with more than 8 elements: ", file_name)
  }
}

# replace the file names to make them more unified 
for (file in final_files) {
  file_name <- basename(file)
  new_file_name <- gsub("alcoholic_hepatitis", "alcoholic-hepatitis", file_name)
  new_file_name <- gsub("lungcancer_copd", "lungcancer_COPD", new_file_name)
  new_file_name <- gsub("adacent", "adjacent", new_file_name)
  new_file_name <- gsub("crc", "CRC", new_file_name)
  new_file_name <- gsub("crc_crc", "CRC_CRC", new_file_name)
  new_file_name <- gsub("CRC_carcinoma", "CRC_CRC", new_file_name)
  new_file_name <- gsub("normal_colon.txt", "normal-colon.txt", new_file_name)
  new_file_name <- gsub("normal_intestine.txt", "normal-intestine.txt", new_file_name)
  new_file_name <- gsub("normal_bronchial_all.txt", "normal-bronchial-all.txt", new_file_name)
  new_file_name <- gsub("Lungcancer", "lungcancer", new_file_name)
  new_file_name <- gsub("normal_cooln.txt", "normal-colon.txt", new_file_name)
  new_file_name <- gsub("normal_N_Cy3.txt", "normal-N-Cy3.txt", new_file_name)
  new_file_name <- gsub("norml_GPL8300.txt", "normal-GPL8300.txt", new_file_name)
  new_file_name <- gsub("normal_ileum.txt", "normal-ileum.txt", new_file_name)
  new_file_name <- gsub("adjacent_old.txt", "adjacent-old.txt", new_file_name)
  new_file_name <- gsub("adjacent_young.txt", "adjacent-young.txt", new_file_name)
  new_file_name <- gsub("pca_pca_vs_", "PCa_PCa_vs_", new_file_name)
  new_file_name <- gsub("pancreaticcancer", "PDAC", new_file_name)
  new_file_name <- gsub("_invasive_cacinoma_ipmn_derived", "_ipmn", new_file_name)
  new_file_name <- gsub("ipmn", "IPMN", new_file_name)
  new_file_name <- gsub("GSE155907_HCC", "GSE155907_livercancer", new_file_name)
  new_file_name <- gsub("GSE143318_HCC", "GSE143318_livercancer", new_file_name)
  new_file_name <- gsub("GSE28619_HCC", "GSE28619_livercancer", new_file_name)
  new_file_name <- gsub("livercancer_HCC_vs", "livercancer_livercancer_vs", new_file_name)
  new_file_name <- gsub("bladdercancer_tumor_vs", "bladdercancer_bladdercancer_vs", new_file_name)
  new_file_name <- gsub("bladdercancer_Ta_T1_vs_normal.txt", "bladdercancer_bladdercancer_vs_normal.txt", new_file_name)
  new_file_name <- gsub("_GPL", "-GPL", new_file_name)
  new_file_path <- file.path(outdir, new_file_name)

  if (new_file_name == file_name) {
    # already correctly named, nothing to do
  } else if (file.exists(new_file_path) && tolower(new_file_name) != tolower(file_name)) {
    message("File already exists with the new name: ", new_file_name)
  } else {
    # macOS/OneDrive is case-insensitive, so a pure case-change rename
    # (e.g. crc_crc -> CRC_CRC) makes file.exists(new_file_path) return
    # TRUE for the *same* file, causing file.rename to be skipped above.
    # Routing through a temporary name avoids that false collision.
    tmp_path <- file.path(outdir, paste0("__tmp_rename__", new_file_name))
    file.rename(file, tmp_path)
    file.rename(tmp_path, new_file_path)
    message("Renamed: ", file_name, " to ", new_file_name)
  }
}

