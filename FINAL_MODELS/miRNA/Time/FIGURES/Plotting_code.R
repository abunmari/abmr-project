source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/Z_level_boxplot.R")
source("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/FIGURES/corr_plots.R")
library(yaml)
library(qgraph)

yaml_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Time/input.yaml"
colors = c("#EECA59", "#CA4F73")
folder = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Time/FIGURES/"

################ LF PLOTS (corr networks and boxplots) ################
# Z11, Z18 , Z8 

l0 = "Early"
l1 = "Late"
group = "Time"

# Z11 plots----------------------------------------------------------------
Z11_plot = Z_level_bp(yaml = yaml_path, Z_num = 11, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z11_boxplot.pdf"), Z11_plot, height = 10, width = 10, scale=0.5)

Z11_corr = plotShape(yaml_path, 11, 0.4, 0.7,  label0 = l0, label1 = l1,
                    save = paste0(folder, "Z11_corr_network_thresh0.4"))

Z11_corr_2 = plotShape(yaml_path, 11, 0.4, 0.25,  label0 = l0, label1 = l1,
                      save = paste0(folder, "Z11_corr_network_thresh0.25"))
# Z18 plots------------------------------------------------------
Z18_plot = Z_level_bp(yaml = yaml_path, Z_num = 18, thresh = c(-2,2),label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z18_boxplot.pdf"), Z18_plot, height = 10, width = 10, scale=0.5)

Z18_corr = plotShape(yaml_path, 18, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z18_corr_network_thresh0.4"))

Z18_corr_2 = plotShape(yaml_path, 18, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z18_corr_network_thresh0.25"))

# Z8 plots------------------------------------------------------
Z8_plot = Z_level_bp(yaml = yaml_path, Z_num = 8, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z8_boxplot.pdf"), Z8_plot, height = 10, width = 10, scale=0.5)

Z8_corr = plotShape(yaml_path, 8, 0.4, 0.7,  label0 = l0, label1 = l1,
                     save = paste0(folder, "Z8_corr_network_thresh0.4"))

Z8_corr_2 = plotShape(yaml_path, 8, 0.4, 0.25,  label0 = l0, label1 = l1,
                       save = paste0(folder, "Z8_corr_network_thresh0.25"))

################ SLIDE BOXPLOT ######################################
source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/slide_bp.R")

SLIDE_plot <- slide_bp("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/miRNA/Time/SLIDE/", TRUE, folder = folder)
ggsave(file = paste0(folder, "SLIDE_performance_boxplot.pdf"), SLIDE_plot, height = 10, width = 10, scale=0.5)


################ CROSS VALIDATION PLOT ################################

source("/ix/djishnu/Marisa/ABomics/cross_preds/cp_SLIDE.R")
canada_mirna_time_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/mirna_31/mirna_x.csv"
canada_mirna_time_y = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/mirna_31/mirna_time_y.csv"

pdf(paste0(folder, "CV_cohort_specific.pdf"))
mirna_time <- load_data_cp(yaml_path = yaml_path, 
                        val_x =canada_mirna_time_x ,
                        val_y =canada_mirna_time_y,
                        interactions = FALSE, scale = TRUE)
dev.off()


