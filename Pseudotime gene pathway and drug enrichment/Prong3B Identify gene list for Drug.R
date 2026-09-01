# This script aiming to identify the best cutoff for genelist selection. Below steps performed for each individual cancer (prostate, colorectal, and lung cancer)
# 1. identify genelist with different cutoff
# 2. Run hallmark pathway and drug enrichment

library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)
library(readxl)
library(dplyr)
library(randomcoloR)
library(factoextra)

# Setup and parameters####
rm(list = ls())
# setwd('/Users/yelin.zhao/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Delade dokument - GRP_MDT_CPMed/EARLY T/Prong3B_Yelin')
inputdir <- 'PseudoTime results from Jonathan 20260715'
outdir <- "Output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
Cancer = 'CRC'#'CRC','prostate_cancer','LUAD'
fc_cutoff = 1.5

# Helpers ####
merge_cancer_data <- function(cancer, inputdir) {
  PT_output <- readRDS(paste0(inputdir, '/', cancer, '_pseudotime_significant_adj_p_val_status_OK.rds')) %>% 
    mutate(gene = rownames(.))
  DE_output <- readRDS(paste0(inputdir, '/', cancer, '_DE_results_all.rds')) %>% 
    mutate(gene = rownames(.))
  Merge_df <- PT_output %>% left_join(DE_output, by = "gene")
  rownames(Merge_df) <- Merge_df$gene
  return(Merge_df)
}

#1 Load data ####
df <- merge_cancer_data(cancer = Cancer, inputdir = inputdir)
dim(df)

#2 Select signatures####
## remove PT and DE non-significant_genes
df_filtered = df %>% filter(status == 'OK' & q_value < 0.05 & padj < 0.05 )
dim(df_filtered)

write.csv(df_filtered$gene, file = paste0(outdir, '/', Cancer, '_filtered_sig_genes_noFCcutoff_noMoransI.csv'), row.names = FALSE)

plot(df_filtered$morans_I, df_filtered$log2FoldChange)

##2.1 remove low FC genes ####
df_stage1 <- df_filtered[abs(df_filtered$log2FoldChange) > fc_cutoff, ]

##2.2 Select genes based on Maximum Perpendicular Distance ####
morans_i_sorted <- sort(df_stage1$morans_I, decreasing = TRUE)
n_points <- length(morans_i_sorted)

# scale
x_scaled <- (1:n_points) / n_points
y_scaled <- (morans_i_sorted - min(morans_i_sorted)) / (max(morans_i_sorted) - min(morans_i_sorted))

x_start <- x_scaled[1]
y_start <- y_scaled[1]
x_end <- x_scaled[n_points]
y_end <- y_scaled[n_points]

# calculate perpendicular distance 
# Distance = |(x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)| / sqrt((x2-x1)^2 + (y2-y1)^2)
numerator <- abs((x_end - x_start) * (y_start - y_scaled) - (x_start - x_scaled) * (y_end - y_start))
denominator <- sqrt((x_end - x_start)^2 + (y_end - y_start)^2)
distances <- numerator / denominator

# get the knee index and Moran's I cutoff
knee_index <- which.max(distances)
optimal_morans_cutoff <- morans_i_sorted[knee_index]

print(paste("Maximum Perpendicular Distance Moran's I cutoff:", round(optimal_morans_cutoff, 4)))

final_shortlist <- df_stage1[df_stage1$morans_I >= optimal_morans_cutoff, ]
print(paste("final_shortlist count:", nrow(final_shortlist)))
final_gene_list = final_shortlist$gene


### Eblow plot ####
png(filename = paste0(outdir, '/', Cancer, '_fc',fc_cutoff,'_elbow_plot.png'), width = 800, height = 600)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(1:n_points, morans_i_sorted, 
     type = "l", lwd = 3, col = "darkgray",
     xlab = "Rank of Genes based on Moran'I", 
     ylab = "Moran's I",
     main = "Elbow Plot: Maximum Perpendicular Distance",
     cex.main = 1.2, cex.lab = 1.2)

points(knee_index, optimal_morans_cutoff, col = "red", pch = 19, cex = 1.8)
abline(v = knee_index, col = "red", lty = 2, lwd = 1.5)
abline(h = optimal_morans_cutoff, col = "blue", lty = 2, lwd = 1.5)

legend("topright", 
       legend = c("Ranked Moran's I Decay", 
                  paste0("Optimal Cutoff (I = ", round(optimal_morans_cutoff, 4), ")")),
       col = c("darkgray", "red"), 
       lwd = c(3, NA), pch = c(NA, 19), lty = c(1, NA), bty = "n")
dev.off()

### Dot plot####
png(filename = paste0(outdir, '/', Cancer, '_fc',fc_cutoff,'_dot.png'), width = 800, height = 600)
plot(df_filtered$morans_I, df_filtered$log2FoldChange, 
     pch = 1, col = "gray80", 
     xlab = "df_filtered$morans_I", 
     ylab = "df_filtered$log2FoldChange",
     main = "Data-Driven Feature Extraction",
     cex.main = 1.2, cex.lab = 1.2)
abline(h = fc_cutoff, col = "darkgreen", lty = 2, lwd = 2)
abline(h = -fc_cutoff, col = "darkgreen", lty = 2, lwd = 2)
abline(v = optimal_morans_cutoff, col = "blue", lty = 2, lwd = 2)
points(final_shortlist$morans_I, final_shortlist$log2FoldChange, 
       pch = 19, col = "red", cex = 1) # selected genes
legend("topright", 
       legend = c("Filtered Shortlist", 
                  paste0("|log2FC| Cutoff = ", fc_cutoff),
                  paste0("Moran's I Cutoff = ", round(optimal_morans_cutoff, 4))),
       col = c( "red", "darkgreen", "blue"),
       pch = c(1, 19, NA, NA), 
       lty = c(NA, NA, 2, 2), 
       lwd = c(NA, NA, 2, 2), bty = "n")
dev.off()

#######################################################
# save final gene list into csv file
write.csv(final_gene_list, file = paste0(outdir, '/', Cancer, '_final_gene_list_fc',fc_cutoff,'.csv'), row.names = FALSE)
 
print(paste("final_shortlist count:", nrow(final_shortlist)))

# write content log.txt 
log_path <- paste0(outdir, "/", Cancer, "_log.txt")
con <- file(log_path, open = "a")
writeLines(c(
  paste("Cancer type:", Cancer),
  paste0("Significant gene list count:", dim(df_filtered)[1]),
  paste("Final gene list count add FC cutoff and elbow cutoff:", length(final_gene_list)),
  paste("Moran's I cutoff:", round(optimal_morans_cutoff, 4)),
  paste("Fold change cutoff:", fc_cutoff),
  ""  # blank line between runs
), con = con) #con = paste0(outdir, '/', Cancer, '_log.txt')
close(con)


#3. Hallmark pathway enrichment ####
# Use Prong3B Enrichment.py to run the enrichment for "final_gene_list"

