
############################################################ 
# plotShape function
############################################################ 
plotShape <- function(yaml_path, lf, threshold, repulsion, spec = NULL, label0 = '0', label1= '1' , colors = c("#EECA59", "#CA4F73"), save = NULL) {
  input <- yaml::yaml.load_file(yaml_path)
  x <- read.csv(input$x_path, row.names = 1, header = TRUE)
  
  sigGenes <- readRDS(paste0(input$out_path, "plotSigGenes_data.RDS"))
  
  print(unique(sigGenes$lf_num))
  
  # Get specific lf 
 
  lf_gene_df <- sigGenes[sigGenes$lf_num == lf,]
  # Unique df
  lf_gene_df <-  lf_gene_df %>%
    distinct(names, .keep_all = TRUE)
  
  if("plain" %in% lf_gene_df$is_marginal){
    title = paste0("Z", lf, " (Interaction)")
  }else{  title = paste0("Z", lf, " (Standalone)")
  }
  
  # Unique lf gene names & color
  lf_genes <- lf_gene_df$names
  lf_colors <- lf_gene_df$color
  print(lf_colors)
  color_labels = recode(lf_colors, 'Blue'= label0, 'Red' = label1, 'white'= "None")
  color_labels = factor(color_labels)
  col_pal = c(label0 = colors[1], label1 = colors[2], "None" = 'lightgray')
  sub <- x[colnames(x) %in% lf_genes]
  corr_matrix <- cor(scale(sub), method = 'spearman')

  q <- qgraph(
    corr_matrix,
    layout = "spring",  
    groups = color_labels,
    color = col_pal,
    posCol="#639354",
    negCol=  "#AF8165",
    rescale = T,
    labels = lf_genes,
    label.scale.equal = TRUE,
    fade = T,
    alpha = 0.05,
    border_width = 0.5,
    threshold = threshold,
    repulsion = repulsion,
    label.cex = 1.2, # scalar on label size
    label.color = 'black', # string on label colors
    border.color = '#3b3b3b',
    vsize = 10,
    legend=TRUE,   
    legend.mode="groups",
    label.font = 1,
    title = title
  ) 
  
  if(!is.null(save)){
    q <- qgraph(
      corr_matrix,
      layout = "spring",  
      groups = color_labels,
      color = col_pal,
      posCol="#639354",
      negCol=  "#AF8165",
      rescale = T,
      labels = lf_genes,
      label.scale.equal = TRUE,
      fade = T,
      alpha = 0.05,
      border_width = 0.5,
      threshold = threshold,
      repulsion = repulsion,
      label.cex = 1.2, # scalar on label size
      label.color = 'black', # string on label colors
      #label.prop = 0.9, # proportion of the width of the node that the label scales
      border.color = '#3b3b3b',
      vsize = 9,
      legend=TRUE,   
      legend.mode="groups",
      label.font = 1,
      title = title, 
      filetype = "pdf",
      filename = save
    ) 
  }
  
  # ggsave(file = paste0(out_path,"delta",input$delta,"lambda",input$lambda,"spec", spec, "_boxplot.png"), lambda_boxplot)
  return(q)
}


############################################################ 
# plotShape function
############################################################ 
plotShape_1 <- function(yaml_path, lf, threshold, repulsion, spec = NULL, label0 = '0', label1= '1' , colors = c( "#CA4F73"), save = NULL) {
  input <- yaml::yaml.load_file(yaml_path)
  x <- read.csv(input$x_path, row.names = 1, header = TRUE)
  
  sigGenes <- readRDS(paste0(input$out_path, "plotSigGenes_data.RDS"))
  
  print(unique(sigGenes$lf_num))
  
  # Get specific lf 
  lf_gene_df <- sigGenes[sigGenes$lf_num == lf,]
  # Unique df
  lf_gene_df <-  lf_gene_df %>%
    distinct(names, .keep_all = TRUE)
  
  if("plain" %in% lf_gene_df$is_marginal){
    title = paste0("Z", lf, " (Interaction)")
  }else{  title = paste0("Z", lf, " (Standalone)")
  }
  
  # Unique lf gene names & color
  lf_genes <- lf_gene_df$names
  lf_colors <- lf_gene_df$color
  
  color_labels = recode(lf_colors, 'Blue'= label0, 'Red' = label1)
  color_labels = factor(color_labels)
  sub <- x[colnames(x) %in% lf_genes]
  corr_matrix <- cor(scale(sub), method = 'spearman')
  
  q <- qgraph(
    corr_matrix,
    layout = "spring",  
    groups = color_labels,
    color = colors,
    posCol="#639354",
    negCol=  "#AF8165",
    rescale = T,
    labels = lf_genes,
    label.scale.equal = TRUE,
    fade = T,
    alpha = 0.05,
    border_width = 0.5,
    threshold = threshold,
    repulsion = repulsion,
    label.cex = 1.2, # scalar on label size
    label.color = 'black', # string on label colors
    #label.prop = 0.9, # proportion of the width of the node that the label scales
    border.color = '#3b3b3b',
    vsize = 9,
    legend=TRUE,   
    legend.mode="groups",
    label.font = 1,
    title = title
  ) 
  
  if(!is.null(save)){
    q <- qgraph(
      corr_matrix,
      layout = "spring",  
      groups = color_labels,
      color = col_pal,
      posCol="#639354",
      negCol=  "#AF8165",
      rescale = T,
      labels = lf_genes,
      label.scale.equal = TRUE,
      fade = T,
      alpha = 0.05,
      border_width = 0.5,
      threshold = threshold,
      repulsion = repulsion,
      label.cex = 1.2, # scalar on label size
      label.color = 'black', # string on label colors
      #label.prop = 0.9, # proportion of the width of the node that the label scales
      border.color = '#3b3b3b',
      vsize = 9,
      legend=TRUE,   
      legend.mode="groups",
      label.font = 1,
      title = title, 
      filetype = "pdf",
      filename = save
    ) 
  }
  
  # ggsave(file = paste0(out_path,"delta",input$delta,"lambda",input$lambda,"spec", spec, "_boxplot.png"), lambda_boxplot)
  return(q)
}


############################################################ 
# plotShape function
############################################################ 
plotShape2 <- function(yaml_path, lf, threshold, repulsion, spec = NULL, label0 = '0', label1= '1' , colors = c("#EECA59", "#CA4F73"), save = NULL) {
  input <- yaml::yaml.load_file(yaml_path)
  x <- read.csv(input$x_path, row.names = 1, header = TRUE)
  
  sigGenes <- readRDS(paste0(input$out_path, "plotSigGenes_data.RDS"))
  
  print(unique(sigGenes$lf_num))
  
  # Get specific lf 
  lf_gene_df <- sigGenes[sigGenes$lf_num == lf,]
  # Unique df
  lf_gene_df <-  lf_gene_df %>%
    distinct(names, .keep_all = TRUE)
  
  if("plain" %in% lf_gene_df$is_marginal){
    title = paste0("Z", lf, " (Interaction)")
  }else{  title = paste0("Z", lf, " (Standalone)")
  }
  
  # Unique lf gene names & color
  lf_genes <- lf_gene_df$names
  lf_colors <- lf_gene_df$color
  
  color_labels = recode(lf_colors, 'Blue'= label0, .default = label1)
  color_labels = factor(color_labels)
  
  sub <- x[colnames(x) %in% lf_genes]
  corr_matrix <- cor(scale(sub), method = 'spearman')
  colnames(corr_matrix)<- remove_ending_number_pattern(colnames(corr_matrix))
  
  q <- qgraph(
    corr_matrix,
    layout = "spring",  
    groups = color_labels,
    color = colors,
    posCol="#639354",
    negCol=  "#AF8165",
    rescale = T,
    labels = lf_genes,
    label.scale.equal = TRUE,
    fade = T,
    alpha = 0.05,
    border_width = 0.5,
    threshold = threshold,
    repulsion = repulsion,
    label.cex = 1.5, # scalar on label size
    label.color = 'black', # string on label colors
    #label.prop = 0.9, # proportion of the width of the node that the label scales
    border.color = '#3b3b3b',
    vsize = 11,
    legend=TRUE,   
    legend.mode="groups",
    label.font = 1
   # title = title
  ) 
  
  if(!is.null(save)){
    q <- qgraph(
      corr_matrix,
      layout = "spring",  
      groups = color_labels,
      color = colors,
      posCol="#639354",
      negCol=  "#AF8165",
      rescale = T,
      labels = lf_genes,
      label.scale.equal = TRUE,
      fade = T,
      alpha = 0.05,
      border_width = 0.5,
      threshold = threshold,
      repulsion = repulsion,
      label.cex = 1.4, # scalar on label size
      label.color = 'black', # string on label colors
      #label.prop = 0.9, # proportion of the width of the node that the label scales
      border.color = '#3b3b3b',
      vsize = 14,
      legend=TRUE,   
      legend.mode="groups",
      label.font = 1,
      #title = title, 
      filetype = "pdf",
      filename = save
    ) 
  }
  
  # ggsave(file = paste0(out_path,"delta",input$delta,"lambda",input$lambda,"spec", spec, "_boxplot.png"), lambda_boxplot)
  return(q)
}
