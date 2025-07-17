source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/Z_level_boxplot.R")
source("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Late/FIGURES/corr_plots.R")
library(yaml)
library(qgraph)

yaml_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/input.yaml"
colors = c("#EECA59", "#CA4F73")
folder = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/FIGURES/"

################ LF PLOTS (corr networks and boxplots) ################
# z1 and z12
group = "Status"
l0 = "DSA+AbMR-"
l1 = "DSA+AbMR+"

# Z1 plots----------------------------------------------------------------

    # Z level boxplot
    Z1_plot = Z_level_bp(yaml = yaml_path, Z_num = 1, label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
    ggsave(file = paste0(folder, "Z1_boxplot.pdf"), Z1_plot, height = 10, width = 10, scale=0.5)

    # Correlation Networks
    Z1_corr = plotShape(yaml_path, 1, 0.4, 0.7, label0 = l0, label1 = l1, save = paste0(folder, "Z1_corr_network_thresh0.4"))
    Z1_corr_2 = plotShape(yaml_path, 1, 0.25, 0.8, label0 = l0, label1 = l1, save = paste0(folder, "Z1_corr_network_thresh0.25"))
# Z12 plots------------------------------------------------------

    # Z level boxplot
    Z12_plot = Z_level_bp(yaml = yaml_path, Z_num = 12,label0 = l0, label1 = l1 , colors = colors, label_order = c(l0,l1), group = group)
    ggsave(file = paste0(folder, "Z12_boxplot.pdf"), Z12_plot, height = 10, width = 10, scale=0.5)

    # Correlation Networks
    Z12_corr = plotShape(yaml_path, 12, 0.4, 0.5, label0 = l0, label1 = l1, save = paste0(folder, "Z12_corr_network_thresh0.4"))
    Z12_corr_2 = plotShape(yaml_path, 12, 0.25, 0.8, label0 = l0, label1 = l1, save = paste0(folder, "Z12_corr_network_thresh0.25"))

################ SLIDE BOXPLOT ######################################
source("/ix/djishnu/Marisa/analysis_scripts/SLIDE_plots/slide_bp.R")

    SLIDE_plot <- slide_bp("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/SLIDE/", TRUE, folder = folder)
    
# cv results in this folder
  #  ~/JishnuLab/ABomics_cleaned/FINAL_MODELS/Abomics/Early/SLIDE
      folder = '~/JishnuLab/ABomics_cleaned/FINAL_MODELS/Abomics/Early/FIGURES/'
    
      lambda_boxplot = ggpubr::ggboxplot(data = perRes, x = "method", y = "auc", palette = "aaas",
                                         fill = "method" ) +
        ggpubr::stat_compare_means(label = "p.signif")
    
    ggplot2::ggsave(plot = lambda_boxplot, filename = paste0(folder, "SLIDECV_boxplot.pdf"), height = 6, width = 7)
    #saveRDS(perRes,file=paste0(slide_input$out_path,"SLIDECV_boxplot_data.rds"))
    

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
 
    
