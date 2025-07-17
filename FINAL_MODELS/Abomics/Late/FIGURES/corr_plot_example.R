
library(yaml)
library(qgraph)
library(dplyr)

plotCorrNet <- function(
    yaml_path, lf, threshold, repulsion, 
    spec = NULL, label0 = '0', label1 = '1', 
    colors = c("#EECA59", "#CA4F73"), save = NULL
) {
  # Load input from YAML file
  input <- yaml::yaml.load_file(yaml_path)
  x <- read.csv(input$x_path, row.names = 1, header = TRUE)
  
  # Load significant genes data
  sig_genes <- readRDS(file.path(input$out_path, "plotSigGenes_data.RDS"))
  
  # Display unique latent factors
  print(unique(sig_genes$lf_num))
  
  # Filter data for the specified latent factor
  lf_gene_df <- sig_genes %>% 
    filter(lf_num == lf) %>% 
    distinct(names, .keep_all = TRUE)
  
  # Determine the plot title
  title <- if ("plain" %in% lf_gene_df$is_marginal) {
    paste0("Z", lf, " (Interaction)")
  } else {
    paste0("Z", lf, " (Standalone)")
  }
  
  # Extract unique gene names and colors
  lf_genes <- lf_gene_df$names
  lf_colors <- lf_gene_df$color
  print(lf_colors)
  
  # Define all possible labels and their corresponding colors
  color_mapping <- c('Blue' = label0, 'Red' = label1, 'white' = "None")
  all_color_labels <- c(label0, label1, "None")
  color_palette <- c(label0 = colors[1], label1 = colors[2], "None" = 'lightgray')
  
  # Map colors to labels and ensure all labels are present
  color_labels <- factor(recode(lf_colors, !!!color_mapping), levels = all_color_labels)
  
  # Subset data for selected genes
  sub_data <- x %>% select(any_of(lf_genes))
  corr_matrix <- cor(scale(sub_data), method = 'spearman')
  
  # Define common graph parameters
  common_qgraph_params <- list(
    layout = "spring",
    groups = color_labels,
    color = color_palette,
    posCol = "#639354",
    negCol = "#AF8165",
    rescale = TRUE,
    labels = lf_genes,
    label.scale.equal = TRUE,
    fade = TRUE,
    alpha = 0.05,
    border_width = 0.5,
    threshold = threshold,
    repulsion = repulsion,
    label.cex = 1.2,
    label.color = 'black',
    border.color = '#3b3b3b',
    vsize = 10,
    legend = TRUE,
    legend.mode = "groups",
    label.font = 1,
    title = title
  )
  
  # Create and display the graph
  q <- do.call(qgraph, c(list(corr_matrix), common_qgraph_params))
  
  # Save the graph if a path is specified
  if (!is.null(save)) {
    common_qgraph_params$vsize <- 9
    common_qgraph_params$filetype <- "pdf"
    common_qgraph_params$filename <- save
    q <- do.call(qgraph, c(list(corr_matrix), common_qgraph_params))
  }
  
  return(q)
}


yaml_path = "FINAL_MODELS/Abomics/Late/input.yaml"
folder = "FINAL_MODELS/Abomics/Late/FIGURES/"
colors = c("#EECA59", "#CA4F73")

################ LF PLOTS (corr networks and boxplots) ################################
# Plots for each LF (Z1, Z2, Z7, Z10) ---------
l0 = "DSA+AbMR-" # Y = 0 label
l1 = "DSA+AbMR+" # Y = 1 Label
 
# plotting for lf Z7 ----------------------------------------------------------------
plotCorrNet(yaml_path, 7, 0.5, 0.7, label0 = l0, label1 = l1)


