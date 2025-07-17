library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(gghalves)
library(tidyr)

setwd("~/JishnuLab/ABomics_7_10/pitt_cohort/")

# AB 36 patients metadata,  X matrix and Y matrix
meta_36 <- readRDS("meta_36.rds")
x = read.csv("null_runs/ab_time/x.csv",row.names=1)

# get patient status (R/NR)
status_y = read.csv("null_runs/ab_status/y.csv", row.names = 1)
status_y$Y[status_y$Y==1] <-"DSA+AbMR+"
status_y$Y[status_y$Y==0] <-"DSA+AbMR-"
colors = c( "#CA4F73", "#EECA59")

# get ig values & log norm 
x_ig <- x %>% select(starts_with("IgG_C"))
x_log = log(x_ig)

# make df for plotting
df = data.frame("Y" = status_y$Y, x_log)
desired_order <- c("DSA+AbMR+", "DSA+AbMR-") 

# Reshape the dataframe from wide to long format
df_long <- df %>%
  pivot_longer(
    cols = -Y,  # Select all columns except 'Y'
    names_to = "Feature",  # Create a new column called 'Feature' for feature names
    values_to = "Value"  # Create a new column called 'Value' for the values
  )

df_long$Y <-  factor(df_long$Y, levels = desired_order)

##############################################################################
# Plot with half violins, p-values, and rotated x-axis labels
##############################################################################

p <- ggplot(df_long, aes(x = Feature, y = Value, fill = Y, color = Y)) +
  geom_half_violin(
    aes(split = Y), 
    alpha = 0.5,
    position = "identity",
    side = "l", 
    trim = T
  ) +
  geom_half_violin(
    aes(split = Y),
    alpha = 0.5,
    position = "identity",
    side = "r", 
    trim = T
  ) +
  labs(title = "IgG Titers", x = "", y = "Measured Titer (log)") +
  theme_classic() +
  scale_fill_manual(values = c("#CA4F73", "#eebc59")) +
  scale_color_manual(values = c("#CA4F73", "#eebc59"), guide = "none") +
  stat_compare_means(
    aes(label = ..p.signif..), 
    method = "wilcox.test", 
    label.x.npc = 'middle', 
    label.y = 10, 
    
    #  label.y = max(df_long$Value),
    size = 3
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(p)

# Adjust the y-axis limits using coord_cartesian() ----------------------
# p_truncated <- p + coord_cartesian(ylim = c(7, 10))  # Set limits
# print(p_truncated)
# ggsave(file = paste0("~/JishnuLab/ABomics_7_10/pitt_cohort/0_manuscript_analyses/Figures/","IgG_violinplot.pdf"), p_truncated , height = 4, width = 12)

