source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/Z_level_boxplot.R")
source("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/FIGURES/corr_plots.R")
library(yaml)
library(qgraph)
# AB REJECTORS MODEL 
# /ix/djishnu/Marisa/ABomics/pitt_cohort/null_runs/ab_r/thresh_0.2/D0.1_1_0.1
yaml_path = "/ix/djishnu/Marisa/ABomics/pitt_cohort/null_runs/ab_r/thresh_0.2/D0.1_1_0.1/yaml_params.yaml"
colors = c("#EECA59", "#CA4F73")
folder = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Additional/Rejectors/"

################ LF PLOTS (corr networks and boxplots) ################
# z1 and z12
group = "Time"
l0 = "Early"
l1 = "Late"

remove_ending_number_pattern <- function(s) {
  # Remove the pattern ".[numbers]." at the end of the string
  sub("\\.[0-9]+\\.$", "", s)
}
remove_prefix <- function(string) {
  return(sub("^[^_]*_", "", string))
}
# Test the function with your examples



# Z2 plots----------------------------------------------------------------

# Z level boxplot
Z2_plot = Z_level_bp(yaml = yaml_path, Z_num = 2, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z2_boxplot.pdf"), Z2_plot, height = 10, width = 10, scale=0.5)

# Correlation Networks
Z2_corr = plotShape2(yaml_path, 2, 0.4, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z2_corr_network_thresh0.4"))

# Z10plots------------------------------------------------------

# Z level boxplot
Z10_plot = Z_level_bp(yaml = yaml_path, Z_num = 10,label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z10_boxplot.pdf"), Z10_plot, height = 10, width = 10, scale=0.5)

# Correlation Networks
Z10_corr = plotShape2(yaml_path, 10, 0.4, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z10_corr_network_thresh0.4"))
#Z12_corr_2 = plotShape(yaml_path, 12, 0.25, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z12_corr_network_thresh0.25"))

# Z11 plots------------------------------------------------------

# Z level boxplot
Z11_plot = Z_level_bp(yaml = yaml_path, Z_num = 11,label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
ggsave(file = paste0(folder, "Z11_boxplot.pdf"), Z11_plot, height = 10, width = 10, scale=0.5)

# Correlation Networks
Z11_corr = plotShape2(yaml_path, 11, 0.4, 0.7, label0 = l0, label1 = l1,colors = c("#CA4F73","#CA4F73"),
 save = paste0(folder, "Z11_corr_network_thresh0.4"))
#Z12_corr_2 = plotShape(yaml_path, 12, 0.25, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z12_corr_network_thresh0.25"))

################ SLIDE BOXPLOT ######################################
source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/slide_bp.R")

SLIDE_plot <- slide_bp("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/", TRUE, folder = folder)

################ CROSS VALIDATION PLOT ################################
source("/ix/djishnu/Marisa/ABomics/cross_preds/cp_SLIDE.R")

canada_early_x = "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_x.csv"
canada_early_y =  "/ix/djishnu/Marisa/ABomics/canada_cohort/Data/ABMR/early_y.csv"

pdf(paste0(folder, "CV_cohort_specific.pdf"))
ab_early<- load_data_cp(yaml_path = yaml_path, 
                        val_x = canada_early_x,
                        val_y =canada_early_y,
                        interactions = FALSE, scale = TRUE)
dev.off()



###############################

early_cp_yaml = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/input.yaml"
early_x = "/ix/djishnu/Marisa/ABomics/pitt_cohort/N32_runs/ab_early/x.csv"
early_y =  "/ix/djishnu/Marisa/ABomics/pitt_cohort/N32_runs/ab_early/y.csv"

ab_early_32 <- load_data_cp(yaml_path = early_cp_yaml, 
                            val_x = early_x,
                            val_y =early_y,
                            interactions = FALSE, scale = TRUE)



