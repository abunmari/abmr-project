# Load necessary libraries
library(stats)  # For correlation and pnorm
library(cocor)
library(reshape2)
library(dplyr)
library(ggplot2)
library(circlize)

# Define the color mapping using colorRamp2
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Function to perform pairwise Fisher Z-transformation test using cocor --------------
pairwise_cocor <- function(cor1, cor2, n1, n2) {
  n <- min(n1, n2)  # Use the smaller sample size for both tests
  p_values <- matrix(NA, nrow = nrow(cor1), ncol = ncol(cor1))
  dimnames(p_values) <- dimnames(cor1)
  
  for (i in 1:(nrow(cor1) - 1)) {
    for (j in (i + 1):ncol(cor1)) {
      test_result <- cocor.indep.groups(cor1[i, j], cor2[i, j], n1, n2, alternative = "two.sided")
      p_values[i, j] <- test_result@fisher1925$p.value
      p_values[j, i] <- p_values[i, j]  # Mirror the p-value
    }
  }
  
  return(p_values)
}

# identify sig delta correlations using Z test (from cocor) ---------------------
sig_diff_z_test <- function(data1, data2){
  # Calculate Spearman correlation matrices
  cor_matrix1 <- rcorr(as.matrix(data1), type = "spearman")$r
  cor_matrix2 <- rcorr(as.matrix(data2), type= "spearman")$r
  
  # Number of observations in each dataset
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  # Perform pairwise Fisher Z-transformation test
  p_values_matrix <- pairwise_cocor(cor_matrix1, cor_matrix2, n1, n2)

  # Adjust for multiple testing using Benjamini-Hochberg method
  upper_tri_indices <- upper.tri(p_values_matrix) # only use upper triangle
  p_adj <- p.adjust(p_values_matrix[upper_tri_indices], method = "BH")
  adjusted_p_values <- matrix(NA, nrow = nrow(p_values_matrix), ncol = ncol(p_values_matrix))
  adjusted_p_values[upper_tri_indices] <- p_adj
  dimnames(adjusted_p_values) <- dimnames(p_values_matrix)  # add dim names 
  return(list(cor_matrix1, cor_matrix2, p_values_matrix, adjusted_p_values ))
}

# helper function to convert matrix to df --------------------
matrix_to_df<- function(mat) {
    # Convert matrix to dataframe and retain row and column names
    df <- as.data.frame(as.table(mat))
    # Rename the columns for clarity
    colnames(df) <- c("Row", "Column", "Value")
    return(df)
  }
  
# Function to filter matrices based on indices from a specific matrix
filter_results <- function(mat_list, matrix_index, threshold) {
    # Get the matrix to use for filtering
    filter_matrix <- mat_list[[matrix_index]]
    # Find the row and column indices where the condition is met (values > 0 and < p-threshold)
    filter_indices <- which((filter_matrix < threshold & filter_matrix > 0), arr.ind = TRUE)
   
     # Get unique row and column indices that meet the condition
    unique_rows <- unique(filter_indices[, 1])
    unique_cols <- unique(filter_indices[, 2])
  
    # Function to filter all matrices in the list based on the filter indices
    filter_all_matrices <- function(mat_list, unique_rows, unique_cols) {
      lapply(mat_list, function(mat) {
        # Apply the row and column indices to each matrix, keeping row and column names
        filtered_mat <- mat[unique_rows, unique_cols, drop = FALSE]
        # Set values greater than the threshold to NA
        return(filtered_mat)
      })
    }
    # Apply the filtering function to all matrices in the list
    filtered_matrices <- filter_all_matrices(mat_list, unique_rows, unique_cols)
    filtered_matrices[[matrix_index]][filtered_matrices[[matrix_index]] > threshold] <- NA
    return(filtered_matrices)
  }
  
# Function to combine matrices to df
combine_matrices_to_dataframe <- function(cor_matrix, pval_matrix) {
  if(length(cor_matrix)== 0){
    print("No sig values to plot")
    return()
  }
    # Convert correlation matrix to dataframe
    cor_df <- as.data.frame(as.table(cor_matrix))
    colnames(cor_df) <- c("Row", "Column", "Correlation")
    
    # Convert p-value matrix to dataframe
    pval_df <- as.data.frame(as.table(pval_matrix))
    colnames(pval_df) <- c("Row", "Column", "PValue")
    
    # Merge the two dataframes by Row and Column
    combined_df <- merge(cor_df, pval_df, by = c("Row", "Column"))
    return(combined_df)
  }
  
# function to make dotplot
make_dotplot <- function(c_mat, p_mat){
    df <- combine_matrices_to_dataframe(c_mat, p_mat)
    # Plotting the dot plot using ggplot2
    dplot <- ggplot(df, aes(x = Column, y = Row)) +
      geom_point(aes(size = PValue, color = Correlation)) + # Size based on p-value, color based on correlation
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0 ,limits =c(-1, 1)) + 
      scale_size(range = c(2, 10)) +  # Adjust the size range for better visibility
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "right") + # Rotate x-axis labels for better readability
      labs(title = "Dot Plot of Spearman Correlation and Adjusted P-Values",
           x = "",
           y = "",
           color = "Spearman Correlation",
           size = "Adjusted P-Value")  # Add labels for the plot
    dplot
  }
  
### rename matrix by removing the extra number
rename_mat <- function(x) {
    colnames(x)<- remove_ending_number_pattern(colnames(x))
    rownames(x)<- remove_ending_number_pattern(rownames(x))
    return(x)
  }
  
# function for pval stars
  labs.function <- function(x) {
    case_when(
      is.na(x) ~ "",       # Return an empty string for NA values
      x >= 0.05 ~ "",      # No asterisk for p-values >= 0.05
      x < 0.05 & x >= 0.01 ~ "*",    # One asterisk for p-values < 0.05 and >= 0.01
      x < 0.01 & x >= 0.001 ~ "**",  # Two asterisks for p-values < 0.01 and >= 0.001
      x < 0.001 ~ "***"    # Three asterisks for p-values < 0.001
    )
  }
  
# plot corr matrix 1 and 2 heatmap with significant corrs, option to remove non-sig
plot_hm_star <-function(mtx_list, sig = F){
    f_m <- lapply(mtx_list, rename_mat)

    # Define the color mapping using colorRamp2
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    
    corr_mat = f_m[[1]]
    corr_mat2 = f_m[[2]]
    p_mat = f_m[[4]]
    
    probes_c = remove_suffix(colnames(corr_mat))
    probes_r = remove_suffix(rownames(corr_mat))
    
    # Apply the function to each element of the p-value matrix
    asterisk_matrix <- apply(p_mat, MARGIN = c(1, 2), labs.function)
    
    if(sig == T){
      corr_mat[p_mat > 0.05 |is.na(p_mat)] <- NA
      corr_mat2[p_mat > 0.05|is.na(p_mat)] <- NA
    }  
      
      # Define function to draw text annotations on the heatmap, ignoring NAs
      draw_asterisks <- function(j, i, x, y, width, height, fill) {
        if (!is.na(asterisk_matrix[i, j])) {  # Check if the value is not NA
          grid.text(asterisk_matrix[i, j], x, y, gp = gpar(fontsize = 12, col = "black"))
        }
      }
    
    # Create Heatmap with Asterisks
    heatmap <- Heatmap(
      corr_mat,
      col = col_fun,
      name = "Spearman \n correlation",
      row_split = probes_r,
      column_split = probes_c, 
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      na_col = "#ededed",
      rect_gp = gpar(col = "black", lwd = 1),
      cell_fun = draw_asterisks , # Add the custom function for asterisks
      row_title_rot = 0
    ) 
    
    # Create Heatmap with Asterisks
    heatmap2 <- Heatmap(
      corr_mat2,
      col = col_fun,
      name = "Spearman\n correlation",
      row_split = probes_r,
      column_split = probes_c, 
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      na_col = "#ededed",
      rect_gp = gpar(col = "black", lwd = 1),
      cell_fun = draw_asterisks , # Add the custom function for asterisks
      row_title_rot = 0
    )
    
    # Create custom legend for p-value significance levels
    p_value_legend <- Legend(
      labels = c("* (p < 0.05)", "** (p < 0.01)", "*** (p < 0.001)"),
      title = "Adjusted P-value \n Significance",
      legend_gp = gpar(fontsize = 10),
      ncol = 1
    )
    
    combined_heatmaps <- heatmap + heatmap2 
    
    # Combine the heatmap and the custom legend, and draw everything
    draw(combined_heatmaps, ht_gap = unit(5, "mm"), annotation_legend_list = list(p_value_legend))  # Adjust the gap as needed (e.g., 5mm)
  }
 
##################
plot_pheatmap <- function(corr_mat){
  # Define a color palette from blue to red
  # Customize and plot the correlation matrix
  col_labels = remove_suffix(colnames(corr_mat))
  row_labels = remove_suffix(row.names(corr_mat))
  
  # Create a data frame for column and row annotations
  col_annotation <- data.frame(Category = factor(col_labels))
  row_annotation <- data.frame(Category = factor(row_labels))
  
  # Define colors for the annotations
  annotation_colors <- list(Category = setNames(brewer.pal(length(unique(col_labels)), "Set2"), unique(col_labels)))
  
  # Determine where to place gaps
  col_gaps <- which(diff(as.numeric(col_annotation$Category)) != 0)
  row_gaps <- which(diff(as.numeric(row_annotation$Category)) != 0)
  
  # Define color palette and breaks
  # color_palette <-  rev(colorRampPalette(brewer.pal(7, "RdYlBu"),  na.value = "#ededed")(200))
  
  
  # Define color palette and breaks
  color_palette <- rev(colorRampPalette(brewer.pal(7, "RdYlBu"))(200))
  breaks <- seq(-1, 1, length.out = 201)  # Ensure the range covers -1 to 1
  
  # Generate the heatmap with annotations
  p <- pheatmap(corr_mat, 
                color = color_palette,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation", 
                border_color = "white",
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_col = col_annotation,
                annotation_row = row_annotation,
                annotation_colors = annotation_colors,
                gaps_col = col_gaps,
                gaps_row = row_gaps,
                breaks = breaks, # range covers -1 to 1
                #na_color = "#ededed",  # Set NA values to light gray
                 display_numbers = p_mat,   
                number_color = "black"  # Color of the numbers
  )
  return(p)
}

