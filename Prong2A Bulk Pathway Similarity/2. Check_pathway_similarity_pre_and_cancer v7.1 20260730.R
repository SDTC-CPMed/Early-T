#### Because more pairs were added and potentially duplicate files exits, and different data structure 
#### This version will first check these files and generate a clearer pool of pathways files to run.
# run this script after run '1.Clean pathway files from Xinxiu.R'

# This version loop through each cancer type and generate heatmap for each cancer type separately, and add more error handling and summary output
# Add p values for similarity measures in the pairwise comparison section
# v7: Add disease_tissue_major_level (broader disease grouping) for heatmap annotation and for the similarity analysis
# v7.1: Only keep GSE datasets flagged UKBB_CT_sig == "Y" in the GEO tracking spreadsheet

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
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(clusterProfiler)  
library(DOSE)
library(randomcoloR)
library(factoextra)
# remove all variables
rm(list = ls())

inputdir <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/Prong2/Prong2A_Final_clean_version_from_YZ/output/Cleaned IPA results'
outdir <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/Prong2/Prong2A_Final_clean_version_from_YZ/output/PathwaySimilarity'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Master GEO dataset tracking sheet - only GSE datasets flagged UKBB_CT_sig == "Y" here are kept
geo_tracking_file <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Mac/Documents/KI-Projects/DiseaseTrajectory/EarlyT/Prong2/Prong2A_Final_clean_version_from_YZ/Bulk_GEO_final_xinxiu yelin update cleaned_VALIDATED_FILTERED.xlsx'


# Helpers ####
clean_ingenuity_numbers <- function(column) {
  # 1. Replace the "minus" sign (unicode) with a standard hyphen
  column <- gsub("−", "-", column)
  
  # 2. Replace European commas with decimal dots
  column <- gsub(",", ".", column)
  
  # 3. Handle the "×10^01" notation: 
  # This regex looks for '×10^' and replaces it with 'E'
  column <- gsub("×10\\^", "E", column)
  
  # 4. Remove any spaces that might have crept in
  column <- gsub("\\s", "", column)
  
  # 5. Convert to numeric
  return(as.numeric(column))
}

# Mapping of each specific disease_tissue value to a broader disease_tissue_major_level
# category, e.g. Adenoma / Adenoma-HGD / adenoma-or-other / Polyp / sessile-serrated-adenoma-polyp
# all collapse to "Adenoma"; CD / UC / IBD collapse to "CD/UC/IBD"; MSI-CRC / CRC collapse to "CRC".
# NOTE: mapping is based on the disease_tissue string alone (regardless of which CT_towards_cancer
# folder it came from), so e.g. disease_tissue == "PDAC" maps to "PDAC" whether it appears under
# CT_towards_cancer "PDAC" or "PanIN". Add new disease_tissue values here as new datasets come in -
# any value not found in this table keeps its original name (with a warning) so nothing silently drops.
disease_tissue_major_map <- c(
  # CRC precursors / cancer
  "Adenoma"                         = "Adenoma/Polyp",
  "Adenoma-HGD"                     = "Adenoma/Polyp",
  "adenoma-or-other"                = "Adenoma/Polyp",
  "Polyp"                           = "Adenoma/Polyp",
  "sessile-serrated-adenoma-polyp"  = "Adenoma/Polyp",
  "CD"                              = "CD/UC/IBD",
  "UC"                              = "CD/UC/IBD",
  "IBD"                             = "CD/UC/IBD",
  "CRC"                             = "CRC",
  "MSI-CRC"                         = "CRC",

  # Cervical
  "CIN3"                            = "CIN",
  "Cervixcancer"                    = "Cervixcancer",

  # Prostate
  "BPH"                             = "BPH",
  "PCa"                             = "PCa",
  "all-PCa"                         = "PCa",

  # Pancreatic
  "Apancreatitis"                   = "Pancreatitis",
  "Apancreatitis-mild"              = "Pancreatitis",
  "Apancreatitis-severe"            = "Pancreatitis",
  "pancreatitis"                    = "Pancreatitis",
  "T2D"                             = "T2D",
  "IPMN"                            = "IPMN",
  "PDAC"                            = "PDAC",

  # Bladder
  "bladdercancer"                   = "bladdercancer",

  # Heart
  "IHD"                             = "IHD",

  # Liver
  "Cirrhosis"                       = "Cirrhosis",
  "cirrhosis-dysplasia-high"        = "Cirrhosis",
  "cirrhosis-dysplasia-low"         = "Cirrhosis",
  "cirrhosis-simple"                = "Cirrhosis",
  "hcv-cirrhosis"                   = "Cirrhosis",
  "Hepatitis"                       = "Hepatitis",
  "alcoholic-hepatitis"             = "Hepatitis",
  "hepatitisb"                      = "Hepatitis",
  "fibrosis-high"                   = "Fibrosis",
  "HCC"                             = "Livercancer",
  "hcc"                             = "Livercancer",
  "hcc-grade1"                      = "Livercancer",
  "hcc-grade2"                      = "Livercancer",
  "hcc-grade3"                      = "Livercancer",
  "early-hcc"                       = "Livercancer",
  "livercancer-hbv"                 = "Livercancer",
  "livercancer-hcv"                 = "Livercancer",
  "malignant-tumor"                 = "Livercancer",
  "livercancer"                     = "Livercancer",

  # Lung
  "lungcancer"                      = "Lungcancer",
  "lungcancer-stageI"               = "Lungcancer",
  "lungcancer-stageII"              = "Lungcancer",
  "COPD"                            = "COPD",
  "asthma"                          = "Asthma",

  # Single-category cancers (kept as-is)
  "pancreaticcancer"                = "pancreaticcancer",
  "uterinecancer"                   = "uterinecancer"
)

get_disease_major_level <- function(disease_tissue_vec) {
  major <- disease_tissue_major_map[disease_tissue_vec]
  unmapped <- unique(disease_tissue_vec[is.na(major)])
  if (length(unmapped) > 0) {
    warning("disease_tissue value(s) not found in disease_tissue_major_map, keeping original name: ",
            paste(unmapped, collapse = ", "))
  }
  major[is.na(major)] <- disease_tissue_vec[is.na(major)]
  unname(major)
}

#1. prepare pathway files####
pathway_list = list.files(inputdir,pattern = "Pathway_")

##1.0 keep only GSE datasets flagged UKBB_CT_sig == "Y" in the GEO tracking sheet####
# Some GSE ids carry a tissue suffix in the pathway filename (e.g. "GSE194331-blood") that isn't
# present in the tracking sheet's "GEO session" column, so strip anything from the first "-" onward
# before matching. NOTE: filtering is at the GSE (dataset) level, matching how the request was framed -
# a handful of GSEs (e.g. GSE47460, GSE63678, GSE148355, GSE19650) have mixed Y/N across different
# sub-comparisons in the tracking sheet, so a GSE is kept here if ANY of its rows is flagged "Y".
geo_tracking <- read_excel(geo_tracking_file, sheet = "Yelin Merged and cleaned")
ukbb_sig_GSE <- geo_tracking %>%
  filter(toupper(trimws(UKBB_CT_sig)) == "Y") %>%
  pull(`GEO session`) %>%
  trimws() %>%
  unique()

pathway_GSE <- sub("-.*", "", sapply(strsplit(pathway_list,"_"),'[[',4))
keep_file <- pathway_GSE %in% ukbb_sig_GSE

cat("UKBB_CT_sig filter: keeping", sum(keep_file), "of", length(pathway_list), "pathway files\n")
cat("Dropped GSE(s) not flagged UKBB_CT_sig == 'Y' in the tracking sheet:\n")
print(sort(unique(pathway_GSE[!keep_file])))

pathway_list <- pathway_list[keep_file]

fc_threshold = sapply(strsplit(pathway_list,"_"),'[[',3)
dataset = sapply(strsplit(pathway_list,"_"),'[[',4)
CT_towards_cancer = sapply(strsplit(pathway_list,"_"),'[[',5)
disease_tissue = sapply(strsplit(pathway_list,"_"),'[[',6)
disease_tissue_major_level = get_disease_major_level(disease_tissue)
comparison = paste0(sapply(strsplit(pathway_list,"_"),'[[',6), "_", "VS", "_",sapply(strsplit(pathway_list,"_"),'[[',8))
comparison = gsub(".txt","",comparison)
index_combined = paste0(CT_towards_cancer,"_",dataset,"_",comparison,"_",fc_threshold)

total_pathways = matrix(NA, ncol=14) %>% as.data.frame()
colnames(total_pathways) = c("Ingenuity Canonical Pathways","-logp","Ratio","z-score", "Molecules",
                             "comparison", "dataset","CT_towards_cancer","disease_tissue","disease_tissue_major_level","fc_threshold","index_combined","adj.p","-log.adj.p")
for (i in 1:length(pathway_list)){
  print(pathway_list[i])
  dd = read.csv(paste0(inputdir,"/",pathway_list[i]), header=F,sep = "\t")
  # dd = read.csv(paste0(inputdir,"/",'Pathway_DEGs_logFC1_GSE41657_CRC_Adenoma-LGD_vs_normal.txt'), header=F,sep = "\t")
  colnames(dd) = c("Ingenuity Canonical Pathways","-logp","Ratio","z-score", "Molecules" )
  dd = dd[2:dim(dd)[1],1:(dim(dd)[2]-1)]
  dd$comparison = comparison[i]
  dd$dataset = dataset[i]
  dd$CT_towards_cancer = CT_towards_cancer[i]
  dd$disease_tissue = disease_tissue[i]
  dd$disease_tissue_major_level = disease_tissue_major_level[i]
  dd$fc_threshold   = fc_threshold[i]
  dd$index_combined = index_combined[i]
  head(dd)
  
  dd$`-logp`   <- clean_ingenuity_numbers(dd$`-logp`)
  dd$Ratio     <- clean_ingenuity_numbers(dd$Ratio)
  dd$`z-score` <- clean_ingenuity_numbers(dd$`z-score`)
  dd$p <- 10^(-dd$`-logp`)
  
  dd$adj.p = p.adjust(dd$p, method="BH")
  dd$p = NULL
  dd$`-log.adj.p` = -log10(dd$adj.p)
  
  total_pathways=rbind(total_pathways,dd)
}
total_pathways = total_pathways[2:dim(total_pathways)[1],] %>% unique()
head(total_pathways)
dim(total_pathways)
table(total_pathways$index_combined)  
write.table(total_pathways,file = paste0(outdir,"/total_pathways_IPA.txt"),row.names =F,quote = F,sep = "\t")

# get the number of unique pathways per dataset
n_ptw_per_dataset = total_pathways %>%
  group_by(index_combined) %>%
  summarise(n = n_distinct(`Ingenuity Canonical Pathways`)) %>% print(n = 1000)
write.table(n_ptw_per_dataset,file = paste0(outdir,"/number_of_unique_pathways_per_dataset.txt"),row.names =F,quote = F,sep = "\t")

##1.2 filter out non-sig pathways####
# total_pathways = total_pathways[total_pathways$`-logp` > -log10(0.05), ]
total_pathways = total_pathways[total_pathways$`-log.adj.p` > -log10(0.05), ]
class(total_pathways)
total_pathways[total_pathways=="NaN"] <- 0 #make z-score as 0 for pathways that enriched but no z-score
head(total_pathways)
dim(total_pathways)
lapply(total_pathways,class)

# check whether there are rows identical consider `Ingenuity Canonical Pathways` and 'index_combined' columns combined
total_pathways %>%
  group_by(`Ingenuity Canonical Pathways`, index_combined) %>%
  summarise(count = n()) %>%
  filter(count > 1)

##1.3 filter out the same dataset but different FC cutoffs####
# since all files has normal range of pathways, then remove the duplicate file that only differs in the fc threshold, and keep the one with fc0.5
unique_index = unique(total_pathways$index_combined) 
cleaned_index <- data.frame(full_name = unique_index) %>%
  extract(full_name, 
          into = c("base_name", "fc_val"), 
          regex = "(.*)_logFC([0-9.]+)$", 
          remove = FALSE) %>%
  mutate(fc_val = as.numeric(fc_val))
#Group by the base name and keep only the row with the max fc_val
unique_index2 <- cleaned_index %>%
  group_by(base_name) %>%
  filter(fc_val == max(fc_val)) %>%
  ungroup()
# final index to keep
final_index <- unique_index2$full_name
print(final_index)

total_pathways = total_pathways %>% filter(index_combined %in% final_index)

# # filter out total_pathways with less than 100 pathways 
# pathway_counts = total_pathways %>%
#   group_by(index_combined) %>%
#   summarise(n_pathways = n()) %>%
#   filter(n_pathways >= 100)
# 
# total_pathways = total_pathways %>% filter(index_combined %in% pathway_counts$index_combined)
#                           
##1.4 Manual removed dataset####
remove_data = c('PDAC_GSE119794_PDAC_VS_normal_logFC0.25', # less than 1000 pathways
                'PDAC_GSE76894_T2D_VS_normal_logFC0.25', # not enough z-scores
                'PCa_GSE88808_PCa_VS_adjacent-old_logFC0.5', # young and old are merged in 'PCa_GSE88808_PCa_adjacent_logFC0.5'
                'PCa_GSE88808_PCa_VS_adjacent-young_logFC0.5',# young and old are merged in 'PCa_GSE88808_PCa_adjacent_logFC0.5'
                'PCa_GSE60371_PCa_VS_normal_logFC1', # data already present using 'PCa_GSE60371_PCa_normal-N-Cy3_logFC0.25',different platform
                'PCa_GSE60371_PCa_VS_normalNCy3_logFC0.25',# data already present using 'PCa_GSE60371_PCa_normal-N-Cy3_logFC0.25', duplicate name
                'PCa_GSE80609_PCa_VS_BPH_logFC1', # data part of 'PCa_GSE80609_all-PCa_VS_BPH_logFC1'
                'PCa_GSE80609_advanced-PCa_VS_BPH_logFC1' ,# data part of 'PCa_GSE80609_all-PCa_VS_BPH_logFC1'
                'lungcancer_GSE19804_ILD_VS_normal_logFC0.5', # wrong file name, exactly same as 'lungcancer_GSE19804_Lungcancer_adjacent_logFC0.5'
                'lungcancer_GSE41861_asthma_VS_normal-bronchial-all_logFC0.25', # only one dataset and not enough z-scores
                'lungcancer_GSE47460_ILD_VS_normal_logFC0.5', # same data but another platform has similar but slightly more pathways 'lungcancer_GSE47460_ILD_normal-GPL6480_logFC0.5'
                'lungcancer_GSE10667_ILD_VS_normal_logFC0.5',# Remove all ILD
                'lungcancer_GSE110147_ILD_VS_normal_logFC1' , # Remove all ILD
                'lungcancer_GSE19804_ILD_VS_normal_logFC0.5', # Remove all ILD
                'lungcancer_GSE21411_ILD_VS_normal-GPL570_logFC0.25', # Remove all ILD
                'lungcancer_GSE24206_ILD_VS_normal_logFC0.5', # Remove all ILD
                'lungcancer_GSE31934_ILD_VS_normal_logFC0.25', # Remove all ILD
                'lungcancer_GSE32537_ILD_VS_normal_logFC0.5', # Remove all ILD
                'lungcancer_GSE47460_ILD_VS_normal-GPL6480_logFC0.5', # Remove all ILD
                'lungcancer_GSE47460_ILD_VS_normal_logFC0.5', # Remove all ILD
                'lungcancer_GSE31210_Lungcancer-stageI_VS_normal_logFC0.5', # data was included in lungcancer_GSE31210_lungcancer_VS_normal_logFC0.5
                'lungcancer_GSE31210_Lungcancer-stageII_VS_normal_logFC0.5', # data was included in lungcancer_GSE31210_lungcancer_VS_normal_logFC0.5
                'lungcancer_GSE179156_lungcancer_asthma_vs_normal-GPL570.txt', # only one dataset and not enough z-scores
                'livercancer_GSE14323_CirrhosisHCC_VS_normal_logFC0.5', # cirrhosis tissue from the same patient as 'livercancer_GSE14323_livercancer_normal_logFC0.5'
                'livercancer_GSE14520_livercancer_VS_adjacent_logFC0.5', # same data but another platform has similar but slightly more pathways 'livercancer_GSE14520_livercancer_adjacent-GPL571_logFC0.5'
                'CRC_GSE75214_CD_VS_normal-ileum_logFC0.25', # same dataset but compare to normal same from different location, similar results so drop this but keep 'CRC_GSE75214_CD_normal-colon_logFC0.25' 
                'CRC_GSE10714_Adenoma_VS_normal_logFC0.25', # Too few pathways with z-score
                'CRC_GSE41657_Adenoma-LGD_VS_normal_logFC1', # same study, similar to HGD, will create a full adenoma vs. normal instead of two grades
                'CRC_GSE37364_Adenoma-LGD_VS_normal_logFC0.5', # same study, similar to HGD, will create a full adenoma vs. normal instead of two grades
                'bladdercancer_GSE3167_tumor_VS_normal_logFC0.5', #same study, similar results so drop this one but keep 'bladdercancer_GSE3167_bladdercancer_normal_logFC0.5' 
                'livercancer_GSE28619_alcoholic-hepatitis_VS_normal_logFC0.25')
total_pathways = total_pathways %>% filter(index_combined %in% remove_data == FALSE)

write.table(total_pathways,file = paste0(outdir,"/total_pathways_IPA_postQC.txt"),row.names =F,quote = F,sep = "\t")

# ======================================================================
# 2 Run heatmap for each unique cancer type####
# ======================================================================
# Get unique cancer types
CT_towards_unique_cancer <- unique(total_pathways$CT_towards_cancer)
cat("Found", length(CT_towards_unique_cancer), "cancer types:\n")
print(CT_towards_unique_cancer)

# Function to generate heatmap for a specific cancer type
generate_cancer_heatmap <- function(cancer_type, pathway_data, output_dir) {
  
  cat("\n", paste(rep("=", 60), collapse = ""))
  cat("\nProcessing cancer type:", cancer_type, "\n")
  
  # Filter data for this cancer type
  cancer_data <- pathway_data[pathway_data$CT_towards_cancer == cancer_type, ]
  
  # Check if there's enough data
  if (nrow(cancer_data) == 0) {
    cat("No enough data for", cancer_type, "- skipping\n")
    return(NULL)
  }
  
  if (length(unique(cancer_data$index_combined)) < 2) {
    cat("Less than 3 unique index_combined for", cancer_type, "- skipping\n")
    return(NULL)
  }
  
  cat("Number of rows:", nrow(cancer_data), "\n")
  cat("Unique index_combined:", length(unique(cancer_data$index_combined)), "\n")
  
  # Prepare the matrix for the heatmap
  heatmap_data <- cancer_data %>%
    dplyr::select(`Ingenuity Canonical Pathways`, index_combined, `z-score`) %>%
    pivot_wider(names_from = index_combined, values_from = `z-score`) %>%
    as.data.frame()
  
  heatmap_data[is.na(heatmap_data)] <- 0  
  
  # Set row names and convert to numeric matrix
  rownames(heatmap_data) <- heatmap_data$`Ingenuity Canonical Pathways`
  heatmap_data$`Ingenuity Canonical Pathways` <- NULL
  mat <- as.matrix(heatmap_data) 
  
  # Remove rows/columns with all zeros
  mat <- mat[rowSums(abs(mat)) > 0, , drop = FALSE]
  if (ncol(mat) > 0) {
    mat <- mat[, colSums(abs(mat)) > 0, drop = FALSE]
  }
  
  write.csv(mat, file = file.path(output_dir, paste0("heatmap_matrix_", gsub(" ", "_", cancer_type), ".csv")))
  
  # Check if matrix has sufficient dimensions
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    cat("Insufficient data for", cancer_type, "- need at least 2 pathways and 2 samples\n")
    cat("Rows:", nrow(mat), "Columns:", ncol(mat), "\n")
    return(NULL)
  }
  
  cat("Matrix dimensions after filtering:", dim(mat)[1], "x", dim(mat)[2], "\n")
  
  # Find optimal k (testing between 2 and 10 clusters)
  opt_k <- fviz_nbclust(mat, kmeans, method = "silhouette", k.max = 10)
  n_row_clusters <- as.numeric(opt_k$data$clusters[which.max(opt_k$data$y)])
  
  # Get metadata for the columns remaining in our matrix
  anno_df <- pathway_data %>%
    dplyr::select(index_combined, dataset, CT_towards_cancer, comparison, disease_tissue, disease_tissue_major_level) %>%
    distinct() %>%
    filter(index_combined %in% colnames(mat)) %>%
    arrange(match(index_combined, colnames(mat)))

  # Make sure order matches matrix columns
  if (!identical(anno_df$index_combined, colnames(mat))) {
    cat("Warning: Column order mismatch - reordering\n")
    mat <- mat[, anno_df$index_combined]
  }

  # Define colors
  dataset_cols <- setNames(rainbow(length(unique(anno_df$dataset))), unique(anno_df$dataset))
  comp_cols <- setNames(topo.colors(length(unique(anno_df$comparison))), unique(anno_df$comparison))
  disease_tissue_cols <- setNames(
    distinctColorPalette(length(unique(anno_df$disease_tissue))),
    unique(anno_df$disease_tissue)
  )
  disease_tissue_major_cols <- setNames(
    distinctColorPalette(length(unique(anno_df$disease_tissue_major_level))),
    unique(anno_df$disease_tissue_major_level)
  )

  # Column annotation
  col_anno <- HeatmapAnnotation(
    Dataset = anno_df$dataset,
    Comparison = anno_df$comparison,
    Disease_tissue = anno_df$disease_tissue,
    Disease_tissue_major_level = anno_df$disease_tissue_major_level,
    col = list(
      Dataset = dataset_cols,
      Comparison = comp_cols,
      Disease_tissue = disease_tissue_cols,
      Disease_tissue_major_level = disease_tissue_major_cols
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )
  
  # Define color mapping
  max_z <- max(abs(mat))
  
  col_fun <- colorRamp2(c(-max_z, 0, max_z), c("blue", "grey", "red"))
  
  # # Determine number of row clusters (try 3 but reduce if needed)
  # n_row_clusters <- min(4, floor(nrow(mat) / 3))
  # if (n_row_clusters < 1) n_row_clusters <- 1
  
  # Generate heatmap
  tryCatch({
    ht <- Heatmap(mat, 
                  name = "z-score",
                  col = col_fun,
                  na_col = "black",
                  
                  # Clustering
                  cluster_rows = FALSE,
                  cluster_columns = TRUE,
                  clustering_distance_rows = 'spearman',#euclidean
                  clustering_distance_columns = 'spearman',# 
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",#complete,ward.D2,average
                  
                  # Row splitting
                  row_km = n_row_clusters,
                  
                  # Annotations
                  top_annotation = col_anno,
                  
                  # Labels and aesthetics
                  show_row_names = FALSE, 
                  show_column_names = TRUE,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 8),
                  column_title = paste0("Clustered IPA Pathways - CT towards ", cancer_type),
                  
                  # Legend customization
                  heatmap_legend_param = list(
                    title = "z-score",
                    at = c(-max_z, 0, max_z),
                    labels = c("Downregulated", "Not Significant or No z", "Upregulated"),
                    color_bar = "continuous",
                    footer = "Black = Non-Significant"
                  )
    )
    
    ht_withpathwayname <- Heatmap(mat, 
                  name = "z-score",
                  col = col_fun,
                  na_col = "black",
                  
                  # Clustering
                  cluster_rows = FALSE,
                  cluster_columns = TRUE,
                  clustering_distance_rows = 'spearman',
                  clustering_distance_columns = 'spearman',
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  
                  # Row splitting
                  row_km = n_row_clusters,
                  
                  # Annotations
                  top_annotation = col_anno,
                  
                  # Labels and aesthetics
                  show_row_names = TRUE, 
                  show_column_names = TRUE,
                  row_names_gp = gpar(fontsize = 3), 
                  column_names_gp = gpar(fontsize = 8),
                  column_title = paste0("Clustered IPA Pathways - CT towards ", cancer_type),
                  
                  # Legend customization
                  heatmap_legend_param = list(
                    title = "z-score",
                    at = c(-max_z, 0, max_z),
                    labels = c("Downregulated", "Not Significant or No z", "Upregulated"),
                    color_bar = "continuous",
                    footer = "Black = Non-Significant"
                  )
    )
    
    # Save to PDF
    output_file <- file.path(output_dir, paste0("heatmap_pathway_similarity_", 
                                                gsub(" ", "_", cancer_type), ".pdf"))
    pdf(output_file, width = 8, height = 11)
    set.seed(123)
    draw(ht)
    dev.off()
    cat("✓ Heatmap saved to:", output_file, "\n")
    
    output_file2 <- file.path(output_dir, paste0("heatmap_pathway_similarity_", 
                                                gsub(" ", "_", cancer_type), "_withpathwayname.pdf"))
    pdf(output_file2, width = 11, height = 40)
    set.seed(123)
    draw(ht_withpathwayname)
    dev.off()
    cat("✓ Heatmap saved to:", output_file2, "\n")
    
    # Also save as PNG for quick viewing (optional)
    png_file <- file.path(output_dir, paste0("heatmap_pathway_similarity_", 
                                             gsub(" ", "_", cancer_type), ".png"))
    png(png_file, width = 11, height = 11, units = "in", res = 300)
    set.seed(123)
    draw(ht)
    dev.off()
    cat("✓ PNG saved to:", png_file, "\n")
    
    
    
  }, error = function(e) {
    cat("Error generating heatmap for", cancer_type, ":\n")
    cat("  ", e$message, "\n")
    return(NULL)
  })
}

# ======================================================================
## MAIN EXECUTION: Run for each cancer type ####
# ======================================================================

# Create summary file
summary_file <- file.path(outdir, "heatmap_generation_summary.txt")
cat("Heatmap Generation Summary\n", file = summary_file)
cat("==========================\n\n", file = summary_file, append = TRUE)
cat("Date:", date(), "\n", file = summary_file, append = TRUE)
cat("\nCancer types processed:\n", file = summary_file, append = TRUE)

# Loop through each cancer type
success_count <- 0
for (cancer_type in CT_towards_unique_cancer) {
  # generate_cancer_heatmap <- function(cancer_type, pathway_data, output_dir)
  result <- generate_cancer_heatmap(cancer_type, total_pathways, outdir)
  
  if (!is.null(result)) {
    success_count <- success_count + 1
    cat(cancer_type, "- SUCCESS\n", file = summary_file, append = TRUE)
  } else {
    cat(cancer_type, "- FAILED/INSUFFICIENT DATA\n", file = summary_file, append = TRUE)
  }
}

# Add summary statistics
cat("\n", file = summary_file, append = TRUE)
cat("==========================\n", file = summary_file, append = TRUE)
cat("Summary:\n", file = summary_file, append = TRUE)
cat("  Total cancer types:", length(CT_towards_unique_cancer), "\n", file = summary_file, append = TRUE)
cat("  Successful:", success_count, "\n", file = summary_file, append = TRUE)
cat("  Failed:", length(CT_towards_unique_cancer) - success_count, "\n", file = summary_file, append = TRUE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Analysis complete!\n")
cat("Results saved to:", outdir, "\n")
cat("Summary file:", summary_file, "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# ======================================================================
# 3 Obtain similarity matrix with P-values ####
# ======================================================================
# --- 1. Identify all unique disease tissues and their associated cancer types ---
# We need a clean mapping of Disease to CT_towards_cancer for the grouping logic
# NOTE: similarity analysis below groups by disease_tissue_major_level (not the finer-grained
# disease_tissue) so e.g. all Adenoma variants are analyzed together, CD/UC/IBD are analyzed
# together, and MSI-CRC/CRC are analyzed together as "CRC".
disease_to_ct <- total_pathways %>%
  dplyr::select(disease_tissue_major_level, CT_towards_cancer) %>%
  distinct() %>% arrange(CT_towards_cancer)

# Pre-define comparison logic ---
# Define which precursor diseases (Dis1) to compare against which Target Cancer (Dis2)
# Keywords here refer to disease_tissue_major_level values
comparison_targets <- list(
  "CRC" = c("CD/UC/IBD", "Adenoma"),#,"IHD"
  "Cervixcancer" = c("CIN"),
  "PCa" = c("BPH"),
  "PDAC" = c("Pancreatitis", "T2D"),
  "lungcancer" = c("COPD"), #"ILD", "asthma"
  "livercancer" = c("Hepatitis", "Cirrhosis")
)

# --- 2. Pre-calculate Meta-Vectors (Unchanged) ---
get_meta_vector <- function(df, disease_name) {
  k = df %>% filter(disease_tissue_major_level == disease_name) %>% pull(dataset) %>% unique() %>% length()

  df %>%
    filter(disease_tissue_major_level == disease_name) %>%
    group_by(`Ingenuity Canonical Pathways`) %>%
    summarise(meta_z = sum(`z-score`, na.rm = TRUE) / sqrt(k)) %>%
    ungroup()
}

all_diseases <- unique(total_pathways$disease_tissue_major_level)
meta_list <- lapply(all_diseases, function(d) get_meta_vector(total_pathways, d))
names(meta_list) <- all_diseases

# --- 3. Pairwise Comparison Loop  ---
results_list <- list()

# Loop through the target cancers defined in our list
for(target_cancer in names(comparison_targets)) {
  
  # Get the list of precursors we want to compare against this cancer
  precursors <- comparison_targets[[target_cancer]]
  
  # Identify the specific 'target' name used in your data (e.g., 'CRC' might be 'CRC_adjacent')
  # This finds any disease_tissue string that contains the target keyword
  actual_target_names <- all_diseases[grep(target_cancer, all_diseases, ignore.case = TRUE)]
  
  for(d2 in actual_target_names) {
    # Now find which precursors exist in our actual data
    for(pre_keyword in precursors) {
      actual_precursor_names <- all_diseases[grep(pre_keyword, all_diseases, ignore.case = TRUE)]
      
      for(d1 in actual_precursor_names) {
        message("Comparing: ", d1, " vs ", d2)
        # Avoid comparing the same thing (e.g., CRC vs CRC)
        if(d1 == d2) next
        
        m1 <- meta_list[[d1]]; m2 <- meta_list[[d2]]
        if(is.null(m1) | is.null(m2)) next
        
        # Align and merge
        comp_df <- tibble(`Ingenuity Canonical Pathways` = unique(c(m1$`Ingenuity Canonical Pathways`, m2$`Ingenuity Canonical Pathways`))) %>%
          left_join(m1, by = "Ingenuity Canonical Pathways") %>%
          left_join(m2, by = "Ingenuity Canonical Pathways") %>%
          mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
          filter(!(meta_z.x == 0 & meta_z.y == 0))
        
        if(nrow(comp_df) > 3) {
          v1 <- comp_df$meta_z.x
          v2 <- comp_df$meta_z.y
          
          # Calculations
          p_test <- cor.test(v1, v2, method = "pearson")
          s_test <- cor.test(v1, v2, method = "spearman", exact = FALSE)
          
          # Cosine
          denom <- (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
          cos_sim <- if(denom != 0) sum(v1 * v2) / denom else 0
          
          # Cosine Permutation P-value
          null_cos <- replicate(1000, {
            v2_shuf <- sample(v2)
            d_shuf <- (sqrt(sum(v1^2)) * sqrt(sum(v2_shuf^2)))
            if(d_shuf == 0) 0 else sum(v1 * v2_shuf) / d_shuf
          })
          cos_p <- sum(abs(null_cos) >= abs(cos_sim)) / 1000
          
          res <- data.frame(
            Precursor_Dis1 = d1,
            Target_Cancer_Dis2 = d2,
            Target_Group = target_cancer,
            Cosine = cos_sim,
            Cosine_P = cos_p,
            Spearman = s_test$estimate,
            Spearman_P = s_test$p.value,
            Pearson = p_test$estimate,
            Pearson_P = p_test$p.value,
            N_Pathways = nrow(comp_df),
            N_Pathways_with_Z_Dis1 = sum(v1 != 0),          # pathways that has z-score in Dis1
            N_Pathways_with_Z_Dis2 = sum(v2 != 0),          # pathways that has z-score in Dis2
            N_Pathways_in_both = sum(v1 != 0 & v2 != 0), # pathways that has z-score in both diseases
            pct_Pathways_in_both = sum(v1 != 0 & v2 != 0) / nrow(comp_df), # percentage of pathways with no z-score in both diseases
            N_Pathways_same_direction = sum(v1 * v2 > 0),   # number of pathways that have same direction (both up or both down)
            pct_Pathways_same_direction = sum(v1 * v2 > 0) / sum(v1 != 0 & v2 !=  0), # percentage of pathways with same direction among those that have z-scores in both
            stringsAsFactors = FALSE
          )
          results_list[[paste(d1, d2, sep="_")]] <- res
        }
      }
    }
  }
}

final_results <- bind_rows(results_list)

# Add Relationship Labels and adjust for multiple testing
final_results <- final_results %>%
  mutate(
    Spearman_Adj_P = p.adjust(Spearman_P, method = "BH"),
    Pearson_Adj_P = p.adjust(Pearson_P, method = "BH"),
    Cosine_Adj_P = p.adjust(Cosine_P, method = "BH")
    
    # Relationship_Type = case_when(
    #   Disease1 == "BPH" ~ "2. Co-incidental",
    #   Target_Group %in% c("CRC", "Cervixcancer") ~ "1a. Proliferative",
    #   Target_Group %in% c("lungcancer", "livercancer", "PDAC") ~ "1b. Inflammation/damage",
    #   TRUE ~ "Unknown"
    # )
  ) %>% arrange( Target_Group, -Cosine)

write.table(final_results, file = paste0(outdir, "/pairwise_similarity_with_pvalues.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")
View(final_results)


