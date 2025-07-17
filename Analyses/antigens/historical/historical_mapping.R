# Load necessary libraries
library(ggplot2)
library(dplyr)
library(lubridate)  # for date handling
library(readxl)
###########################################################################################################
# Function to split and trim names
split_and_trim <- function(names_str, remove = NULL) {
  names <- unlist(strsplit(names_str, ","))  
  names <- trimws(names)
  names <- names[!is.na(names)]
  if (!is.null(remove)) {
    names <- names[!(names %in% remove)]
  }
  return(names)
}

# Function to create presence/absence matrix
create_presence_absence_matrix <- function(name_list, remove = NULL) {
  # Split and trim names
  split_names <- lapply(name_list, split_and_trim, remove)
  
  # Get all unique names
  all_names <- unique(unlist(split_names))
  
  # Initialize matrix
  presence_matrix <- matrix(0, nrow = length(name_list), ncol = length(all_names))
  rownames(presence_matrix) <- names(name_list)
  colnames(presence_matrix) <- all_names
  
  # Fill the matrix
  for (i in seq_along(name_list)) {
    presence_matrix[i, colnames(presence_matrix) %in% split_names[[i]]] <- 1
  }
  
  return(presence_matrix)
}

########################################################################################

### based on serotype ####
# Group the data by ID and then extract the names into lists
# Apply the function to each row in the data frame and group names by ID
# names_list <- lapply(split(hist$SEROTYPE, hist$PAT_ID), unique)
# names_list <- lapply(split(hist$SEROTYPE, hist$PAT_ID), 
#                      function(names_str) unique(unlist(lapply(names_str, split_and_trim))))

pat_ag <- function(x){
  names_list <- lapply(split(x$SEROTYPE, x$PAT_ID), unique)
  names_list <- lapply(split(x$SEROTYPE, x$PAT_ID), 
                       function(names_str) unique(unlist(lapply(names_str, split_and_trim))))
  
}
########################################################################################

hist <- read_excel("~/JishnuLab/ABomics_cleaned/pitt_cohort/antigens/historical/PITT_Historical_DSA_Full_5_20_trirupa.xlsx",
                   sheet = 1)

within30_tx <- hist[hist$DAYS_TX_TO_RESULT <= 30,]
within30_list <- pat_ag(within30_tx)


within30_samp <- hist[hist$`DAYS_ SAMPLE_TO_RESULT` <= 30,]
within30_samp_list <- pat_ag(within30_samp)
length(within30_samp)


post30_tx <- hist[hist$DAYS_TX_TO_RESULT > 30,]
post30_list <- pat_ag(post30_tx)

within90_hist <- hist[hist$SAMPLE_WITHIN_90_DAYS == "Yes",]
within90_list <- pat_ag(within90_hist)

within90_tx <- hist[hist$DAYS_TX_TO_RESULT < 90 & hist$DAYS_TX_TO_RESULT > -90,]
withintx_list <- pat_ag(within90_tx)
withintx_list <- match_fill_lists(post90_list, withintx_list)

post90 <- hist[hist$SAMPLE_WITHIN_90_DAYS == "No",]
post90_list <- pat_ag(post90)

find_unique_and_shared_per_sample <- function(list1, list2) {
  # Check if both lists have the same length
  if (length(list1) != length(list2)) {
    if(length(list1)>length(list2)){
      names = list1
      fill = list2
    }else{
      names = list2
      fill = list1}
    fill = match_fill_lists(names,fill )
  }
  
  # Use lapply to apply the same function to each pair of samples
  result <- lapply(seq_along(names), function(i) {
    vector1 <- names[[i]]
    vector2 <- fill[[i]]
    
    shared <- intersect(vector1, vector2)
    unique1 <- setdiff(vector1, vector2)
    unique2 <- setdiff(vector2, vector1)
    
    list(l1 = unique1, l2 = unique2, shared = shared)
  })
  
  names(result) <- names(names)
  return(result)
}

result <- find_unique_and_shared_per_sample(within90_list, post90_list)
                                
result2 <- find_unique_and_shared_per_sample(withintx_list,  post90_list)

result3 <- find_unique_and_shared_per_sample(within30_list,  post30_list)


hist_filt <- hist[!is.na(hist$SEROTYPE),]
