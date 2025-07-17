# creating eGFR plots for pitt & toronto cohort
# plot for manuscript
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stats)
library(cowplot)
library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts

# read in pitt metadata
pitt_meta_path = "~/JishnuLab/ABomics_7_10/pitt_cohort/meta_36.rds"
meta_36 <- readRDS("~/JishnuLab/ABomics_7_10/pitt_cohort/meta_36.rds")
clin_feats = c("Tacrolimus.at.sample", "eGFR.at.Bx", "eGFR.at.sample" , "MVI.g.ptc.2."  , "Interstitial.inflam.i.score") 
colors_picked = c("#CA4F73",  "#eebc59")
meta_36$GRP <- recode(meta_36$GRP, "ABMR" = "DSA+AbMR+", "DSA"="DSA+AbMR-")
desired_order <- c("DSA+AbMR+", "DSA+AbMR-") 

############################
# PITSBURGH eGFR
############################
df <- meta_36[,c("GRP", clin_feats[2])]
t1 <- gsub(">60", 61, df[,2])
t1 <- as.numeric(t1)
df[,2]<- t1

# Remove rows with NA values
df_cleaned <- na.omit(df)
# Convert GRP to a factor with the specified levels
df_cleaned$GRP <- factor(df_cleaned$GRP, levels = desired_order)

p1 <- ggplot(df_cleaned, aes(x = GRP, y = eGFR.at.Bx, 
                             color = GRP,
                             fill = as.factor(GRP))) +
  scale_fill_manual(values = colors_picked, guide = "none") +
  scale_color_manual(values = colors_picked, guide = "none") +
  labs(title = "eGFR at Biopsy",
       x = "",
       y = "eGFR",
       color = "Patient Status") +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(  # RANGE plot 
    adjust = .25, ##  .25 bandwidth
    width = .4, ##  .67 bandwidth
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3
  ) + stat_compare_means(method = "wilcox.test", aes(label =  ..p.signif..), label.x = 1.5, label.y = max(df_cleaned$eGFR.at.Bx) * 1.1) + 
  theme_classic() + theme(panel.grid.major = element_blank(),   
                          plot.title = element_text(hjust = 0.5, face = "bold") ,
                          panel.grid.minor = element_blank()) +guides(fill = guide_legend(ncol = 1, title = "Patient Status"))

############################
# TORONTO eGFR plot
############################
canada_meta_path = "~/JishnuLab/ABomics_7_10/canada_cohort/updated_datasets/canada_meta.rds")
c_meta = readRDS(canada_meta_path)
c_meta$GRP <- c_meta$`Study Group`
# "Number of gloms"                                                 
# %Globally sclerosed gloms"   
c_meta$GRP <- recode( c_meta$GRP,"ABMR" = "DSA+AbMR+", "NR"="DSA+AbMR-")

df <- c_meta[,c("GRP", "Number of gloms")]
df$eGFR<- as.numeric(df$`Number of gloms`)

# Remove rows with NA values
df_cleaned <- na.omit(df)
# Convert GRP to a factor with the specified levels
df_cleaned$GRP <- factor(df_cleaned$GRP, levels = desired_order)

p2 <- ggplot(df_cleaned, aes(x = GRP, y = eGFR, 
                             color = GRP,
                             fill = as.factor(GRP))) +
  scale_fill_manual(values = colors_picked, guide = "none") +
  scale_color_manual(values = colors_picked, guide = "none") +
  labs(title = "eGFR at Biopsy",
       x = "",
       y = "eGFR",
       color = "Patient Status") +
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(  # RANGE plot 
    adjust = .25, ##  .25 bandwidth
    width = .4, ##  .67 bandwidth
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)
  ) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3
  ) + stat_compare_means(method = "wilcox.test", aes(label =  ..p.signif..), label.x = 1.5, label.y = max(df_cleaned$eGFR) * 1.1) + 
  theme_classic() + theme(panel.grid.major = element_blank(),   
                          plot.title = element_text(hjust = 0.5, face = "bold") ,
                          panel.grid.minor = element_blank()) +guides(fill = guide_legend(ncol = 1, title = "Patient Status"))

############################
# COMBINE PLOTS
############################

library(patchwork)


# Create the label for the first plot
label1 <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Pittsburgh Cohort", angle = 90, size = 5, hjust = 0.5) + 
  theme_void() +
  theme(plot.background = element_rect(fill = "lightgrey", linewidth = 3, color = "white"))

# Create the label for the second plot
label2 <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Toronto Cohort", angle = 90, size = 5, hjust = 0.5) + 
  theme_void() +
  theme(plot.background = element_rect(fill = "lightgrey", linewidth = 3, color = "white"))

# Combine plots and labels
combined_plot <- label1 + p1 + label2 +p2  +
  plot_layout(ncol = 2, nrow =2, heights = c(1, 1), widths = c(1,10))

# Display the combined plot

print(combined_plot)
#ggsave(file = paste0("Figures/","egfr_pitt_toronto.pdf"), combined_plot , height = 7, width = 7)

