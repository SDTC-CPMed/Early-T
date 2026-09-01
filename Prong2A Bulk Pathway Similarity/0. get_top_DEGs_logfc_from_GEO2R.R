library(readxl)
library(tidyverse)
library(readr)
library(purrr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)  
library(stringr)
# BiocManager::install("hugene10sttranscriptcluster.db")
# library(hugene10sttranscriptcluster.db)

##test
# AnnotationDbi::columns(hgu133plus2.db)
# AnnotationDbi::keytypes(hgu133plus2.db)


# ---- directories ----
input_dir <- '/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Delade dokument - GRP_MDT_CPMed/EARLY T/Datasets/Additional datasets added by yelin/Raw DEG batch3 to xinxiu'

output_dir <- file.path(input_dir, "toplogFC_output")   # or set another folder you like

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List all TSV files with FULL paths
ls <- list.files(path = input_dir, pattern = "\\.tsv$", full.names = TRUE)
print(ls)

# ---- function to process one file ----
DEG_process <- function(one_file) {
  # print(paste("Processing file:", basename(one_file)))
  
  data <- read.delim(one_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  # ---- Gene symbol mapping (unchanged) ----
  if ("Gene.symbol" %in% colnames(data)) { 
    data <- data %>%
      mutate(GeneID = sub("\\s*///.*$", "", Gene.symbol))
  } else if ("Gene.Symbol" %in% colnames(data)) {
    data <- data %>%
      mutate(GeneID = Gene.Symbol)
  } else if ("Symbol" %in% colnames(data)) {
    data <- data %>%
      mutate(GeneID = Symbol)
  } else if ("gene_assignment" %in% colnames(data)) {
    data <- data %>%
      mutate(GeneID = as.character(gene_assignment))
  } else if ("ENTREZ_GENE_ID" %in% colnames(data)) {
    entrez_ids <- data %>%
      pull(ENTREZ_GENE_ID) %>%
      as.character() %>%
      str_trim() %>%
      discard(~ is.na(.x) || .x == "") %>%
      unique()
    mapped <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL",
                     keytype = "ENTREZID", multiVals = "first")
    data <- data %>%
      mutate(ENTREZ_GENE_ID_chr = str_trim(as.character(ENTREZ_GENE_ID)),
             GeneID = unname(mapped[ENTREZ_GENE_ID_chr])) %>%
      dplyr::select(-ENTREZ_GENE_ID_chr)
  } else if ("SPOT_ID" %in% colnames(data) & "GB_ACC" %in% colnames(data) &
             !("Gene.symbol" %in% colnames(data)) & !("gene_assignment" %in% colnames(data)) &
             !("Gene.Symbol" %in% colnames(data))) {
    data <- data %>%
      mutate(clean_ensembl = str_split_i(as.character(SPOT_ID), ";", 1) %>% str_trim(),
             clean_refseq  = str_split_i(as.character(GB_ACC),  ";", 1) %>% str_trim())
    ensembl_keys <- data$clean_ensembl[!is.na(data$clean_ensembl) & data$clean_ensembl != ""] %>% unique()
    refseq_keys  <- data$clean_refseq[!is.na(data$clean_refseq)  & data$clean_refseq  != ""] %>% unique()
    if (any(grepl("^ENST", ensembl_keys))) {
      mapped_ensembl <- mapIds(org.Hs.eg.db, keys = ensembl_keys, column = "SYMBOL",
                               keytype = "ENSEMBLTRANS", multiVals = "first")
    } else {
      mapped_ensembl <- mapIds(org.Hs.eg.db, keys = ensembl_keys, column = "SYMBOL",
                               keytype = "ENSEMBL", multiVals = "first")
    }
    mapped_refseq <- mapIds(org.Hs.eg.db, keys = refseq_keys, column = "SYMBOL",
                            keytype = "ACCNUM", multiVals = "first")
    data <- data %>%
      mutate(sym_ensembl = unname(mapped_ensembl[clean_ensembl]),
             sym_refseq  = unname(mapped_refseq[clean_refseq]),
             GeneID      = coalesce(sym_ensembl, sym_refseq)) %>%
      dplyr::select(-clean_ensembl, -clean_refseq, -sym_ensembl, -sym_refseq)
  } else if ("SPOT_ID" %in% colnames(data) & !("GB_ACC" %in% colnames(data))) {
    gene_ids <- data %>%
      pull(SPOT_ID) %>%
      as.character() %>%
      str_split_i(";", 1) %>%
      str_trim() %>%
      discard(~ is.na(.x) || .x == "") %>%
      unique()
    mapped <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL",
                     keytype = "ENSEMBL", multiVals = "first")
    data <- data %>%
      mutate(GENE_ID_chr = str_trim(as.character(SPOT_ID)),
             GeneID      = unname(mapped[GENE_ID_chr])) %>%
      dplyr::select(-GENE_ID_chr)
  } else {
    message("Skipping file (no gene symbol column): ", basename(one_file))
    return(NA)
  }
  
  # ---- Standardise p-value / logFC column names ----
  if ("padj"           %in% colnames(data)) data <- data %>% mutate(adj.P.Val = padj)
  if ("log2FoldChange" %in% colnames(data)) data <- data %>% mutate(logFC = log2FoldChange)
  
  if (!all(c("adj.P.Val", "logFC") %in% colnames(data))) {
    message("Skipping file (missing adj.P.Val or logFC): ", basename(one_file))
    return(NA)
  }
  
  # ---- Clean + collapse to one row per gene ----
  data2 <- data %>%
    mutate(GeneID = trimws(GeneID)) %>%
    filter(!is.na(GeneID), GeneID != "") %>%
    group_by(GeneID) %>%
    summarize(adj.P.Val = mean(adj.P.Val, na.rm = TRUE),
              logFC     = mean(logFC,     na.rm = TRUE),
              .groups   = "drop")
  
  # Significant genes only (p < 0.05) — pool used for all threshold checks
  sig_data <- data2 %>% filter(adj.P.Val < 0.05)
  
  # ---- Adaptive logFC threshold ----
  # Count DEGs at each candidate threshold
  thresholds   <- c(0.25, 0.5, 1.0)
  n_per_thresh <- sapply(thresholds, function(thr) sum(abs(sig_data$logFC) > thr))
  names(n_per_thresh) <- as.character(thresholds)
  
  # Selection logic: default 0.5; relax to 0.25 if < 100; tighten to 1 if > 5000
  n_default <- n_per_thresh["0.5"]
  
  chosen_threshold <- dplyr::case_when(
    n_default < 1000  ~ 0.25,
    n_default > 5000 ~ 1.0,
    TRUE             ~ 0.5
  )
  
  # Final filtered set using the chosen threshold
  filtered_data <- sig_data %>%
    filter(abs(logFC) > chosen_threshold) %>%
    arrange(adj.P.Val)
  
  n_chosen <- nrow(filtered_data)
  cat(sprintf("%-55s  thresholds [0.25=%d | 0.5=%d | 1.0=%d]  →  chosen=%.2f (%d DEGs)\n",
              basename(one_file),
              n_per_thresh["0.25"], n_per_thresh["0.5"], n_per_thresh["1"],
              chosen_threshold, n_chosen))
  
  # ---- Write DEG file ----
  file_name   <- tools::file_path_sans_ext(basename(one_file))
  output_file <- file.path(output_dir,
                           paste0("DEGs_logFC", chosen_threshold, "_", file_name, ".tsv"))
  if (n_chosen > 50) {write_tsv(filtered_data, output_file)}
  
  # ---- Return summary row ----
  data.frame(
    dataset          = file_name,
    n_DEGs_thr0.25   = n_per_thresh["0.25"],
    n_DEGs_thr0.5    = n_per_thresh["0.5"],
    n_DEGs_thr1.0    = n_per_thresh["1"],
    chosen_threshold = chosen_threshold,
    n_DEGs_chosen    = n_chosen,
    output_file      = output_file,
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
}


# Apply to all files ####
output_files <- lapply(ls, DEG_process)

deg_counts <- dplyr::bind_rows(output_files)
readr::write_tsv(deg_counts, file.path(output_dir, "/DEG_counts_summary.tsv"))


print(output_files)


###############################
