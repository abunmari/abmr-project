
################################################################
# ~~~~~ Remove ending numbers -----
################################################################
remove_ending_number_pattern <- function(s) {
  # Remove the pattern ".[numbers]." at the end of the string
  sub("\\.[0-9]+\\.$", "", s)
}

################################################################
# ~~~~~ Remove prefix -----
################################################################
remove_prefix <- function(string) {
  return(sub("^[^_]*_", "", string))
}

################################################################
# ~~~~~ Remove suffix -----
################################################################
remove_suffix <- function(string) {
  return(sub("_([^_]*)$", "", string))
}

################################################################
# Function to extract letters from a string
################################################################
extract_letters <- function(string) {
  # Use gsub to remove all non-letter characters
  letters_only <- gsub("[^A-Za-z]", "", string)
  return(letters_only)
}

################################
# Function to split and trim names
################################
split_and_trim <- function(names_str, remove = NULL) {
  names <- unlist(strsplit(names_str, ","))  
  names <- trimws(names)
  names <- names[!is.na(names)]
  if (!is.null(remove)) {
    names <- names[!(names %in% remove)]
  }
  return(names)
}
################################################################
# Function to split strings into alphabetical and numerical parts
################################################################
split_alphanumeric <- function(strings) {
  # Extract letters
  letters_part <- gsub("[^A-Za-z]", "", strings)
  
  # Extract numbers
  numbers_part <- as.numeric(gsub("[^0-9]", "", strings))
  
  # Combine into a data frame
  df <- data.frame(original = strings, letters = letters_part, numbers = numbers_part)
  
  # Order by letters and then by numbers
  df <- df[order(df$letters, df$numbers), ]
  
  return(df$original)
}

##########################################################
# Function to get non-enriched antigens (based on negative)
##########################################################
not_enriched <- function(matrix_data) {
  # Check if each element of the column is non-negative
  non_negative_cols <- apply(matrix_data, 2, function(col) all(col >= 0))
  # Subset the matrix to keep only non-negative columns
  filtered_matrix <- matrix_data[, non_negative_cols]
  return(colnames(filtered_matrix))
}

##########################################################
# Functions to highlight enriched ags in matrix
##########################################################
highlight_enriched_probes <- function(matrix, enrich_df) {
  highlighted_matrix <- matrix
  for (i in seq_len(nrow(enrich_df))) {
    probe <- enrich_df$Probe[i]
    antigen <- enrich_df$QueryAg[i]
    if (probe %in% rownames(matrix) && antigen %in% colnames(matrix)) {
      highlighted_matrix[probe, antigen] <- highlighted_matrix[probe, antigen] * -1
    }
  }
  return(highlighted_matrix)
}

highlight_enriched_probes_filt <- function(matrix, enrich_df) {
  highlighted_matrix <- matrix
  highlight_probes <-remove_suffix(rownames(matrix))
  for (i in seq_len(nrow(enrich_df))) {
    probe <- enrich_df$Probe[i]
    antigen <- enrich_df$QueryAg[i]
    if (probe %in% highlight_probes && antigen %in% colnames(matrix)) {
      h_i <- which(highlight_probes %in% probe)
      highlighted_matrix[h_i, antigen] <- highlighted_matrix[h_i, antigen] * -1
    }
  }
  return(highlighted_matrix)
}

##########################################################
# Function to get all enriched matrices
##########################################################
get_enrich_matrices<-function(matrix, ABMR_enrich, NR_enrich){
  matrix <- as.matrix(matrix)
  highlight_ABMR <- highlight_enriched_probes(matrix, ABMR_enrich)
  n1 <- not_enriched(highlight_ABMR)
  
  highlight_NR <- highlight_enriched_probes(matrix, NR_enrich)
  n2 <- not_enriched(highlight_NR)
  
  remove_feats <- intersect(n1, n2)
  keep <- setdiff(colnames(matrix), remove_feats)
  
  # # Final matrices for plotting
  enrich_mtx <- matrix[, keep]
  highlight_ABMR <- highlight_ABMR[, keep]
  highlight_NR <- highlight_NR[, keep]
  return(list("enrich_mtx"=enrich_mtx, 
              "highlight_ABMR" = highlight_ABMR, 
              "highlight_NR" = highlight_NR))
}

##########################################################
# Function to get all enriched matrices - latent factors
##########################################################
get_enrich_matrices_lf <-function(matrix, ABMR_enrich, NR_enrich){
  matrix <- as.matrix(matrix)
  highlight_ABMR <- highlight_enriched_probes_filt(matrix, ABMR_enrich)
  n1 <- not_enriched(highlight_ABMR)
  
  highlight_NR <- highlight_enriched_probes_filt(matrix, NR_enrich)
  n2 <- not_enriched(highlight_NR)
  
  remove_feats <- intersect(n1, n2)
  keep <- setdiff(colnames(matrix), remove_feats)
  
  # # Final matrices for plotting
  enrich_mtx <- matrix[, keep]
  highlight_ABMR <- highlight_ABMR[, keep]
  highlight_NR <- highlight_NR[, keep]
  return(list("enrich_mtx"=enrich_mtx, 
              "highlight_ABMR" = highlight_ABMR, 
              "highlight_NR" = highlight_NR))
}

##############################################################
# Functions for individual corr plots 
##############################################################

### function for pval stars ##################################
labs.function <- function(x) {
  case_when(
    is.na(x) ~ "",  
    x == 0 ~ "",       # Return an empty string for NA values
    # Return an empty string for NA values
    x >= 0.05 ~ "",      # No asterisk for p-values >= 0.05
    x < 0.05 & x >= 0.01 ~ "*",    # One asterisk for p-values < 0.05 and >= 0.01
    x < 0.01 & x >= 0.001 ~ "**",  # Two asterisks for p-values < 0.01 and >= 0.001
    x < 0.001 ~ "***"    # Three asterisks for p-values < 0.001
  )
}

### corr m ptest ################################################
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

### corrplot fx ################################################
AB_corrplot <- function(data, thresh){
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  M <- cor(data, method = "spearman")
  p.mat <- cor.mtest(data, method = "spearman", exact = NULL)
  # Specialized the insignificant value according to the significant level
  plot <- corrplot(M, type="lower", insig = "label_sig", diag= F,
                   p.ma = p.mat, sig.level = thresh,col = col_fun,
                   # cl.lim=c(-1,0,1), col=colorRampPalette(c("blue","white","red")(200)),
                   tl.col="black", tl.srt=45 )
  plot
}

##################################################################################
AB_heatmap<- function(sub, thresh = NULL){
  corr_mat <- cor(scale(sub), method = 'spearman')
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  probes_c = remove_suffix(colnames(corr_mat))
  probes_r = remove_suffix(rownames(corr_mat))
  p_mat <- cor.mtest(corr_mat, method = "spearman", exact = NULL)
  corr_mat[p_mat > thresh |is.na(p_mat)] <- NA
  
  if(is.null(thresh)){
    
    asterisk_matrix <- apply(p_mat, MARGIN = c(1, 2), labs.function)
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
      #cell_fun = draw_asterisks , # Add the custom function for asterisks
      row_title_rot = 0
    ) 
    heatmap
  }else{ heatmap <- Heatmap(
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
    row_title_rot = 0)
  heatmap
  }
}

################################################
# Misc helper functions
################################################
get_pairs <- function(df){
  df_filtered <- df %>%
    mutate(Row = as.character(Row), Column = as.character(Column)) %>%
    filter(Row <= Column)
  return(df_filtered)
}

convert_to_matrix <- function(df, colname){
  matrix <- acast(df, Row ~ Column, value.var = colname)
  return(matrix)
}


# ##########################################################
# # Function to add missing patients to a list 
# ##########################################################
# match_fill_lists <- function(names_list,list ){
#   expected_keys = names(names_list)
#   # Add missing entries with a value of 0
#   for (key in expected_keys) {
#     if (!key %in% names(list)) {
#       list[[key]] <- "0"
#     }
#   }
#   list <- list[expected_keys]
#   return(list)
# }

# ################################################################
# # ~~~~~ Presence/Absence matrix from list of features -----
# ################################################################
# PA_matrix <- function(list_of_lists) {
#   # Identify unique elements
#   unique_elements <- unique(unlist(list_of_lists))
#   # Initialize the presence-absence matrix
#   presence_absence_matrix <- matrix(0, nrow = length(list_of_lists), ncol = length(unique_elements))
#   colnames(presence_absence_matrix) <- unique_elements
#   row.names(presence_absence_matrix)<- names(list_of_lists)
#   # Fill the matrix
#   for (i in 1:length(list_of_lists)) {
#     for (element in list_of_lists[[i]]) {
#       presence_absence_matrix[i, element] <- 1
#     }
#   }
#   # Return the presence-absence matrix
#   return(presence_absence_matrix)
# }
# 
# 
# ################################################################
# # Function to create presence/absence matrix
# ################################################################
# create_presence_absence_matrix <- function(name_list, remove = NULL) {
#   # Split and trim names
#   split_names <- lapply(name_list, split_and_trim, remove)
#   
#   # Get all unique names
#   all_names <- unique(unlist(split_names))
#   
#   # Initialize matrix
#   presence_matrix <- matrix(0, nrow = length(name_list), ncol = length(all_names))
#   rownames(presence_matrix) <- names(name_list)
#   colnames(presence_matrix) <- all_names
#   
#   # Fill the matrix
#   for (i in seq_along(name_list)) {
#     presence_matrix[i, colnames(presence_matrix) %in% split_names[[i]]] <- 1
#   }
#   
#   return(presence_matrix)
# }
