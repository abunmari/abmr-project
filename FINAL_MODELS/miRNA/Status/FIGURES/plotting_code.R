source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/Z_level_boxplot.R")
source("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/FIGURES/corr_plots.R")
library(yaml)
library(qgraph)

yaml_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Status/input.yaml"
colors = c("#EECA59", "#CA4F73")
folder = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Status/FIGURES/"

################ LF PLOTS (corr networks and boxplots) ################
# 103 111  61  69 135


l0 = "DSA+AbMR-"
l1 = "DSA+AbMR+"
group = "Status"

# Z111 plots----------------------------------------------------------------
Z111_plot = Z_level_bp(yaml = yaml_path, Z_num = 111, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z111_boxplot.pdf"), Z111_plot, height = 10, width = 10, scale=0.5)

Z111_corr = plotShape(yaml_path, 111, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z111_corr_network_thresh0.4"))

Z111_corr_2 = plotShape(yaml_path, 111, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z111_corr_network_thresh0.25"))


# Z103 plots------------------------------------------------------
Z103_plot = Z_level_bp(yaml = yaml_path, Z_num = 103, thresh = c(-2,2),label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z103_boxplot.pdf"), Z103_plot, height = 10, width = 10, scale=0.5)

Z103_corr = plotShape(yaml_path, 103, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z103_corr_network_thresh0.4"))

Z103_corr_2 = plotShape(yaml_path, 103, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z103_corr_network_thresh0.25"))


# Z61 plots------------------------------------------------------
Z61_plot = Z_level_bp(yaml = yaml_path, Z_num = 61, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z61_boxplot.pdf"), Z61_plot, height = 10, width = 10, scale=0.5)

Z61_corr = plotShape(yaml_path, 61, 0.4, 0.7,  label0 = l0, label1 = l1,
                    save = paste0(folder, "Z61_corr_network_thresh0.4"))

Z61_corr_2 = plotShape(yaml_path, 61, 0.4, 0.25,  label0 = l0, label1 = l1,
                      save = paste0(folder, "Z61_corr_network_thresh0.25"))

# Z69 plots------------------------------------------------------
Z69_plot = Z_level_bp(yaml = yaml_path, Z_num = 69, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z69_boxplot.pdf"), Z69_plot, height = 10, width = 10, scale=0.5)

Z69_corr = plotShape(yaml_path, 69, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z69_corr_network_thresh0.4"))

Z69_corr_2 = plotShape(yaml_path, 69, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z69_corr_network_thresh0.25"))


# Z135 plots------------------------------------------------------
Z135_plot = Z_level_bp(yaml = yaml_path, Z_num = 135, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z135_boxplot.pdf"), Z135_plot, height = 10, width = 10, scale=0.5)

Z135_corr = plotShape(yaml_path, 135, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z135_corr_network_thresh0.4"))

Z135_corr_2 = plotShape(yaml_path, 135, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z135_corr_network_thresh0.25"))

################ SLIDE BOXPLOT ######################################
source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/slide_bp.R")

SLIDE_plot <- slide_bp("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Time/SLIDE/", TRUE, folder = folder)
ggsave(file = paste0(folder, "SLIDE_performance_boxplot.pdf"), SLIDE_plot, height = 10, width = 10, scale=0.5)


################ CROSS Pred PLOT ################################

source("/ix/djishnu/Marisa/ABomics/cross_preds/cp_SLIDE.R")
canada_mirna_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/mirna_31/mirna_x.csv"
canada_mirna_y = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/mirna_31/mirna_group_y.csv" 


mirna_status <- load_data_cp(yaml_path = yaml_path, 
                           val_x =canada_mirna_x ,
                           val_y =canada_mirna_y,
                           interactions = TRUE, scale = TRUE)

ggsave(file = paste0(folder, "CV_plot.pdf"), mirna_status, height = 10, width = 10, scale=0.5)

