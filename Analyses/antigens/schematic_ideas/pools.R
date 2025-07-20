# Code to help visualize the minipools
setwd('~/JishnuLab/ABomics_cleaned')
source("Analyses/antigens/plotting_codes/Plotting_helper_fxns.R")

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


# 19 pools
ag_mapping_list <- readRDS("~/JishnuLab/ABomics_7_10/pitt_cohort/antigens/ag_mapping_list.rds")

# viewing the 19 pools as a matrix 
ag_mapping_matrix <- as.matrix(readRDS("~/JishnuLab/ABomics_7_10/pitt_cohort/antigens/ag_mapping_matrix.rds"))


filt_matrix <- ag_mapping_matrix

f_mat <- function(x){
  cols_to_keep <- grepl("^(A|DP)", colnames(x))
  x[,cols_to_keep, drop = FALSE]
}

filt_matrix <- f_mat(filt_matrix)
antigen_types = extract_letters(colnames(filt_matrix))

t <-Heatmap(filt_matrix, 
            col = c("white", "#3A69AE"),
            name = "counts",
            # row_split = probe_label,
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            cluster_columns = FALSE,
            cluster_rows= FALSE,
            column_split = antigen_types,
            row_names_side = "left",
            
            column_gap = unit(5, "mm"),        
            rect_gp = gpar(col = "black", lwd = 1),
            row_names_gp = gpar(fontsize = 11),
            column_names_gp = gpar(fontsize = 11)
)

folder = "~/JishnuLab/ABomics_cleaned/Analyses/antigens/schematic_ideas/"
name = "probe_blue.pdf"
pdf(paste0(folder,name), width = 7, height = 4)  # Specify file name and dimensions
draw(t)  # Plot the heatmap to the PDF
dev.off()  # Close the PDF device
# # # 


full_ag_matrix = readRDS("~/JishnuLab/ABomics_7_10/pitt_cohort/antigens/full_feat_ag_matrix.rds")

t2 <-Heatmap(full_ag_matrix, 
            col = c("white", "grey"),
            name = "counts",
            # row_split = probe_label,
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            cluster_columns = FALSE,
            cluster_rows= FALSE,
            column_split = antigen_types,
            row_names_side = "left",
            
            column_gap = unit(3, "mm"),        
            rect_gp = gpar(col = "black", lwd = 1),
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10)
)