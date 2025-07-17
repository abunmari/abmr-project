
library(EssReg)
library(ggplot2)
library(yaml)
library(dplyr)
library(pROC)

# slide_res$SLIDE_res$interaction_vars
# [1] "Z1.Z5"  "Z1.Z11" "Z12.Z3"
## predZ ################################################################################################################
predZ <- function(x, er_res) {
  A_hat <- er_res$A
  C_hat <- er_res$C
  Gamma_hat <- er_res$Gamma
  Gamma_hat <- ifelse(Gamma_hat == 0, 1e-10, Gamma_hat)
  Gamma_hat_inv <- diag(Gamma_hat ** (-1))
  G_hat <- crossprod(A_hat, Gamma_hat_inv) %*% A_hat + solve(C_hat)
  Z_hat <- x %*% Gamma_hat_inv %*% A_hat %*% MASS::ginv(G_hat)
  return (Z_hat)
}

### Pairwise interaxtions #########################################################
pairwiseInteractions <- function(index_list, mat) {
  num_cols <- ncol(mat)
  index_combinations <- expand.grid(seq_len(num_cols),index_list,stringsAsFactors=F)
  temp <- index_combinations$Var1
  index_combinations$Var1 <- index_combinations$Var2
  index_combinations$Var2 <- temp
  
  col_names <- paste0(colnames(mat)[as.numeric(index_combinations[, 1])], ".", colnames(mat)[as.numeric(index_combinations[, 2])])
  interaction_mat <- mat[, as.numeric(index_combinations[, 1])] * mat[, as.numeric(index_combinations[, 2])]
  if(is.null(dim(interaction_mat))){
    interaction_mat <- matrix(interaction_mat,nrow=1)}
  
  colnames(interaction_mat) <- col_names
  return(list(interaction=as.data.frame(interaction_mat)))
}

### sig LF df #########################################################
getLFdf<- function(SLIDE_res, z_matrix, interactions = TRUE){
  # marginal data
  sigK <- SLIDE_res$marginal_vals
  sigK <- c(sigK)
  sigMargData <- z_matrix[,sigK]
  
  if(!interactions){
    return(as.data.frame(sigMargData))
  }
  
  # interaction
  sigIn <- SLIDE_res$SLIDE_res$interaction_vars
  if(!is.null(sigIn)){
    IntData <- pairwiseInteractions(sigK,z_matrix)
    sigIntData <- data.frame(IntData$interaction[ ,sigIn])
    names(sigIntData)<- sigIn
    
    final_mtx <- cbind(as.data.frame(sigMargData), sigIntData)
    return(final_mtx)
  }
  return(as.data.frame(sigMargData))
}

# cross prediction for SLIDE #########################################################
cross_prediction_SLIDE <- function(er_results, train_x, train_y,  val_x, val_y, slide_res, interactions, scale, save = NULL, lqbels = NULL) {
  
  train_in_val <- colnames(train_x)[which(colnames(train_x) %in% colnames(val_x))]
  # Get column names in train_x that are not in val_x
  train_not_in_val <- colnames(train_x)[which(!(colnames(train_x) %in% colnames(val_x)))]


  if (length(train_in_val) != length(colnames(train_x) )){
    print("Should re-run ER & SLIDE with matched feature set ")
    # Create a data frame with 0 values for the new columns
    new_cols <- data.frame(matrix(0, nrow = nrow(val_x), ncol = length(train_not_in_val)))
    colnames(new_cols) <- train_not_in_val
    val_x <- cbind(val_x, new_cols)
  }
  
  # if validation dataset has aditional features, remove them 
  val_not_in_train <- colnames(val_x)[which(!colnames(val_x) %in% colnames(train_x))]
  val_x <- val_x[, !colnames(val_x) %in% val_not_in_train]
  
  val_x <- as.matrix(val_x)
  train_x <- as.matrix(train_x)
  
  if(scale){
  val_x <- scale(as.matrix(val_x))
  train_x <- scale(as.matrix(train_x))
  }
  
  train_z <- predZ(train_x , er_results)
  colnames(train_z) <- paste0("Z", c(1:ncol(train_z)))
  
  val_z <- predZ(val_x, er_results)
  colnames(val_z) <- paste0("Z", c(1:ncol(val_z)))
  
  # get LF dataframe for link function
  train_input = getLFdf(slide_res, train_z, interactions)
  val_input = getLFdf(slide_res, val_z, interactions)
  
  print("LFs used for model")
  print(colnames(train_input))
  # link function 
  lin_reg <- stats::glm(train_y ~ ., data = train_input, family = "gaussian")
  
  # predicting on train data
  train_y_pred   <- predict(lin_reg, train_input ,type = 'response')
  train_auc <- pROC::auc(response=as.matrix(train_y), predictor=as.matrix(train_y_pred),quite=T)
 
  # predicting on validation data
  val_y_pred   <- predict(lin_reg, val_input ,type = 'response')
  val_auc <- pROC::auc(response=as.matrix(val_y), predictor=as.matrix(val_y_pred),quite=T)
 
  # roc for plotting
  roc_score_train = roc(train_y, train_y_pred)
  roc_score_val = roc(val_y, val_y_pred)
  
  # Plot the ROC curves with enhanced aesthetics
  cp_plot <- plot.roc(roc_score_train, 
                      main = "ROC Curve",
                      colorize = TRUE, 
                      col = "#1687A7",  
                      lwd = 3, 
                      print.auc = TRUE, 
                      print.auc.x = 0.5, 
                      print.auc.y = 0.9,
                      grid = FALSE,  # Add grid lines for better readability
                      grid.col = "white")  # Light gray grid lines for subtlety
  
    # Add the validation ROC curve
    plot.roc(roc_score_val, 
             add = TRUE, 
             colorize = TRUE, 
             col = "#DD0A35",
             lwd = 3,  # Match line width with the train ROC curve
             print.auc = TRUE, 
             print.auc.x = 0.7, 
             print.auc.y = 0.5)
    
    if(is.null(labels)){
      legend =  c("Train", "Validation")
    }else{
      train_label = paste0("Train:", labels[1])
      val_label = paste0("Val: ", labels[2])
      legend=c(train_label, val_label)
    }
    
    # Customize the legend
    legend(x = "bottomright", 
           box.lwd = 1, 
           box.col = "darkgray",  # Dark gray box for a cleaner look
           bg = "white",  # White background for the legend box
           title = "ROC Curves", 
           legend = legend, 
           fill = c( "#1687A7", "#DD0A35"),  # Match legend colors with plot lines
           border = "white")  # No border around the fill colors in the legend
    
  # Display the plot
  cp_plot
  return(list(plot = cp_plot,
              train_auc = round(train_auc, 2),
              val_auc = round(val_auc, 2), 
              model = lin_reg,
              lfs =  colnames(train_input)
              ))
} 

### load data for crossprediction
# marg_z = c (1, 2, 3)
# int_z = c("Z1.Z5", "Z1.Z11")
load_data_cp <- function(yaml_path, val_x, val_y, interactions = TRUE, scale = FALSE, new = NULL, marg_z = NULL, int_z = NULL, labels = NULL){
  
  input <- read_yaml(yaml_path)

  if(!is.null(new)){  # input$out_path <- new
    er_path = list.files(input$out_path, full.names = T,  pattern = "AllLatent")
    er_results = readRDS(er_path)
    slide_path= list.files(input$out_path,recursive = T,  full.names = T,  pattern = "SLIDE_LF")
    slide_res = readRDS(slide_path)
  }else{
    er_path = list.files(input$out_path, full.names = T,  pattern = "final_")
    er_results = readRDS(er_path)
    slide_res = readRDS(list.files(input$out_path,recursive = T, full.names = T, pattern = "slide_res")[1])
  }
  train_x <- as.matrix(read.csv(input$x_path, row.names=1))
  train_y <- as.matrix(read.csv(input$y_path, row.names=1))
  val_x <-as.matrix(read.csv(val_x, row.names=1))
  val_y <- as.matrix(read.csv(val_y, row.names=1))
  

  if(!is.null(marg_z)){
    slide_res$marginal_vals <- marg_z
  }
  if(!is.null(int_z)){
    slide_res$SLIDE_res$interaction_vars <- int_z
  }
  cp <- cross_prediction_SLIDE(er_results, train_x, train_y, val_x, val_y, slide_res, interactions, scale, labels)
  
  return(cp)
}
