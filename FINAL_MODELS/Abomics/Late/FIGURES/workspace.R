library(yaml)
library(qgraph)
setwd("~/JishnuLab/ABomics_cleaned/")
#source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/Z_level_boxplot.R")
source("plotting_codes/SLIDE_corr_networks.R")

yaml_path = "FINAL_MODELS/Abomics/Late/input.yaml"
folder = "FINAL_MODELS/Abomics/Late/FIGURES/"
colors = c("#EECA59", "#CA4F73")
other_colors = c("#DDD9A8")

################ LF PLOTS (corr networks and boxplots) ################################
# Plots for each LF (Z1, Z2, Z7, Z10) ---------
group = "Status"
l0 = "DSA+AbMR-"
l1 = "DSA+AbMR+"

# Z1 ----------------------------------------------------------------

    # Z level boxplot
    Z1_plot = Z_level_bp(yaml = yaml_path, Z_num = 1,thresh = c(-0.75,0), group = group, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1))
    ggsave(file = paste0(folder, "Z1_boxplot.pdf"), Z1_plot, height = 10, width = 10, scale=0.5)
   
    # Correlation Networks
    Z1_corr = plotShape(yaml_path, 1, 0.4, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z1_corr_network"))


# Z2 ----------------------------------------------------------------
    
    # Z level boxplot
    Z2_plot = Z_level_bp(yaml = yaml_path, Z_num = 2, group = group,label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1))
    ggsave(file = paste0(folder, "Z2_boxplot.pdf"), Z2_plot, height = 10, width = 10, scale=0.5)
    
    # Correlation networks
    Z2_corr = plotShape(yaml_path, 2, 0.4, 0.7, label0 = l0, label1 = l1,
                         save = paste0(folder, "Z2_corr_network"))

# Z7 ----------------------------------------------------------------
    
    # Z level boxplot
    Z7_plot = Z_level_bp(yaml = yaml_path, Z_num = 7, group = group, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1))
    ggsave(file = paste0(folder, "Z7_boxplot.pdf"), Z7_plot, height = 10, width = 10, scale=0.5)
    
    # Correlation network
    Z7_corr = plotShape(yaml_path, 7, 0.4, 0.7, label0 = l0, label1 = l1, colors = c("#CA4F73", "lightgray" ),
                        save = paste0(folder, "Z7_corr_network"))

# Z10 ----------------------------------------------------------------

    # Z level boxplot    
    Z10_plot = Z_level_bp(yaml = yaml_path, Z_num = 10, group = group,label0 = l0, label1 = l1 , 
                          colors = colors, label_order = c(l0,l1), thresh = c(-1.2,1.2))
    ggsave(file = paste0(folder, "Z10_boxplot_thresh1.2.pdf"), Z10_plot, height = 10, width = 10, scale=0.5)
    
    # Correlation Network
    Z10_corr = plotShape(yaml_path, 10, 0.4, 0.7, label0 = l0, label1 = l1,
                         save = paste0(folder, "Z10_corr_network"))


################ SLIDE BOXPLOT ######################################
source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/slide_bp.R")

SLIDE_plot <- slide_bp("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/SLIDE/", TRUE, folder = folder)


################ CROSS VALIDATION PLOT ################################

source("/ix/djishnu/Marisa/ABomics/cross_preds/cp_SLIDE.R")
canada_late_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_x.csv"
canada_late_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/late_y.csv"

# DEFINING Y USING THE COHORT SPECIFIC MEDIAN 
pdf(paste0(folder, "CV_cohort_specific.pdf"))
ab_late<- load_data_cp(yaml_path = yaml_path, 
                        val_x = canada_late_x,
                        val_y =canada_late_y,
                        interactions = FALSE, scale = TRUE)

dev.off()



# USING ONLY THE PITT MEDIAN

canada_late_new_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/new_late/AB_late_med_x.csv"
canada_late_new_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/new_late/AB_late_med_y.csv"


pdf( paste0(folder, "CV_pitt_median.pdf"))
ab_late_median <- load_data_cp(yaml_path = yaml_path, 
                                   val_x = canada_late_new_x,
                                   val_y =canada_late_new_y,
                                   interactions = FALSE, scale = TRUE)

dev.off()

#ggsave(file = paste0(folder, "CV_pitt_median.pdf"), ab_late_median, height = 10, width = 10, scale=0.5)




