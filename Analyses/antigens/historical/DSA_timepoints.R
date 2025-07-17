# Load necessary libraries
setwd("~/JishnuLab/ABomics_cleaned/")

library(ggridges)
library(ggplot2)
library(dplyr)
library(lubridate)  # for date handling
library(readxl)
library(dplyr)
source("~/JishnuLab/ABomics_cleaned/Analyses/antigens/plotting_codes/Plotting_helper_fxns.R")
source("Analyses/correlation/Helper_fx.R")
source("Analyses/load_data.R")

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

# get historical data and patient data
hist<- readRDS("~/JishnuLab/ABomics_cleaned/Analyses/antigens/historical/pitt_historical_10_2.rds")
# early and late patients
patient_early_ag_data <- readRDS("~/JishnuLab/ABomics_cleaned/Analyses/antigens/ag_data/patient_early_ag_data.rds")
patient_late_ag_data <- readRDS("~/JishnuLab/ABomics_cleaned/Analyses/antigens/ag_data/patient_late_ag_data.rds")

######################################################
# get late hist only
######################################################
late_id = c(patient_late_ag_data$ID_abmr, patient_late_ag_data$ID_dsa)
hist_late <- hist[hist$PAT_ID %in% late_id,]

# if a serotype is stored as a list, will expand this list into additional rows
late_df_sep<- hist_late %>%
  separate_rows(SEROTYPE, sep = ",") %>%
  drop_na(SEROTYPE)
late_df_sep$SEROTYPE <- trimws(late_df_sep$SEROTYPE)

# get data from within a timeframe
#  patients 501 , 40, 451  have a very late date
pat_cutoff = as.numeric(late_df_sep$DAYS_TX_TO_RESULT) <= 1000
late_df_separated<- late_df_sep[pat_cutoff,]
length(unique(late_df_separated$PAT_ID))

# Remove duplicates and create unique serotypes per day
df_unique <- late_df_separated %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)

# get antigen data
df_unique$type <- extract_letters(df_unique$SEROTYPE)

#####################################
# line plots over time 
#####################################
# Create line plot of serotypes over time in late cohort -------------------------------
time_late <- ggplot(df_unique, aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n (Late group)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 30 days
  
# Create line plot of serotypes over time in late cohort for only REJECTOR patients -------------------------------
df_separated_ABMR <- late_df_separated[late_df_separated$Status == "ABMR",] %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)
df_separated_ABMR$type <- extract_letters(df_separated_ABMR$SEROTYPE)

time_ABMR_late <- ggplot(df_separated_ABMR , aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n (LATE: REJECTORS)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60days

 # PLOT sep by patients test
  # df_separated_ABMRpat <- late_df_separated[late_df_separated$Status == "ABMR",] %>%
  #   distinct(PAT_ID,DAYS_TX_TO_RESULT, SEROTYPE)
  # df_separated_ABMRpat$type <- extract_letters(df_separated_ABMRpat$SEROTYPE)
  # 
  # time_ABMR_latepat <- ggplot(df_separated_ABMRpat , aes(x =DAYS_TX_TO_RESULT, y = factor(PAT_ID), group =PAT_ID, color = type)) +
  #   geom_line(size = 1) + 
  #  #geom_text(hjust=0, vjust=0) +
  #   geom_point(size = 2) +
  #   labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n (LATE: REJECTORS)") +
  #   geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60days
  # 
# Create line plot of serotypes over time in late cohort for only NON-REJECTOR patients --------------
df_separated_DSA <- late_df_separated[late_df_separated$Status == "DSA",] %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)
df_separated_DSA$type <- extract_letters(df_separated_DSA$SEROTYPE)

time_DSA_late <- ggplot(df_separated_DSA , aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n  (LATE: NR)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60days

#########################################
# first DSA positive date for each patient (Dotplots)
#########################################
# group by patients
df_uniquepat<- late_df_separated %>%
  distinct(PAT_ID, DAYS_TX_TO_RESULT, SEROTYPE, Status)
df_uniquepat$type <- extract_letters(df_uniquepat$SEROTYPE)

# get first (shortest) date
df_min <- df_uniquepat %>%
  group_by(PAT_ID) %>%
  filter(DAYS_TX_TO_RESULT == min(DAYS_TX_TO_RESULT)) %>%
  ungroup()

# Plot the minimum days to result for each patient (Scatter Plot)
# 
firstDSA_late <- ggplot(df_min, aes(x = DAYS_TX_TO_RESULT, y =factor(PAT_ID), color = Status,  label = PAT_ID))+
  geom_point() + 
 geom_text(hjust=0, vjust=0) +
  geom_point(size = 3) +
  labs( title = "Days till first DSA result (LATE)") +
  theme_minimal() +   
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60 days
R_ID = c("351","287", "146", "55")
###################################################
# Ridgeplots: Serotype distribution across days 
###################################################

# Create a ridge plot for all late patients
ggplot(df_unique, aes(x = DAYS_TX_TO_RESULT, y = type, height = ..density..)) +
  geom_density_ridges(stat = "density", scale = 2, fill = "lightblue", color = "blue") +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution Across Days (Late patient group)") +
  theme_minimal()  +   
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60 days

# get patient status column
df_unique2 <- late_df_separated %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE, Status)
df_unique2$type<- extract_letters(df_unique2$SEROTYPE)
df_unique2$Status <- factor(df_unique2$Status)

# NR/R patients sep
ggplot(df_unique2, aes(x = DAYS_TX_TO_RESULT, y = type, height = ..density..)) + 
  geom_density_ridges(stat = "density", scale = 2) + 
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution Across Days (Late patient group)") + 
  theme_minimal() + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~Status)

# NR/R patients overlap
late_ridge <-ggplot(df_unique2, aes(x = DAYS_TX_TO_RESULT,fill=Status, y = type, height = ..density..)) + 
  geom_density_ridges(stat = "density", scale = 2, alpha = 0.5) + 
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution Across Days (Late patient group)") + 
  theme_minimal() + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1) 

plot_grid(firstDSA_late, late_ridge, time_late, time_DSA_late, time_ABMR_late)
plot_grid(time_DSA_late, time_ABMR_late)

## NOW AT EARLY TIME POINT
######################################################
# get EARLY hist only
######################################################
early_id = c(patient_early_ag_data$ID_abmr, patient_early_ag_data$ID_dsa)
hist_early <- hist[hist$PAT_ID %in% early_id,]

# if a serotype is stored as a list, will expand this list into additional rows
early_df_sep <- hist_early %>%
  separate_rows(SEROTYPE, sep = ",") %>%
  drop_na(SEROTYPE)
early_df_sep$SEROTYPE <- trimws(early_df_sep$SEROTYPE)

# get data from within a year (n = 16)
early_df_separated<- early_df_sep[early_df_sep$DAYS_TX_TO_RESULT<=4000,]
# Remove duplicates and create unique serotypes per day
df_unique_early <- early_df_separated %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)
# get antigen data
df_unique_early$type <- extract_letters(df_unique_early$SEROTYPE)

#####################################
# line plots over time 
#####################################
# Create line plot of serotypes over time in early cohort -------------------------------
time_early <- ggplot(df_unique_early, aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n(Early Group)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 30 days

# Create line plot of serotypes over time in early cohort for only REJECTOR patients -------------------------------
early_df_separated_ABMR <- early_df_separated[early_df_separated$Status == "ABMR",] %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)
early_df_separated_ABMR$type <- extract_letters(early_df_separated_ABMR$SEROTYPE)

time_ABMR_early <- ggplot(early_df_separated_ABMR , aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n (Early: REJECTORS)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60days

# Create line plot of serotypes over time in early cohort for only NON-REJECTOR patients --------------
early_df_separated_DSA <- early_df_separated[early_df_separated$Status == "DSA",] %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE)
early_df_separated_DSA$type <- extract_letters(early_df_separated_DSA$SEROTYPE)

time_DSA_early <- ggplot(early_df_separated_DSA , aes(x =DAYS_TX_TO_RESULT, y = SEROTYPE, group = SEROTYPE, color = type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Occurrences Over Time \n (Early: NR)") +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60days


# both plots combined
plot_grid(time_DSA_early, time_ABMR_early)

#########################################
# first DSA positive date for each patient (Dotplots)
#########################################
# group by patients
df_unique_earlypat <- early_df_separated %>%
  distinct(PAT_ID, DAYS_TX_TO_RESULT, SEROTYPE, Status)
df_unique_earlypat$type <- extract_letters(df_unique_earlypat$SEROTYPE)

# get first (shortest) date
df_min_early <- df_unique_earlypat %>%
  group_by(PAT_ID) %>%
  filter(DAYS_TX_TO_RESULT == min(DAYS_TX_TO_RESULT)) %>%
  ungroup()

# Plot the minimum days to result for each patient (Scatter Plot)
firstDSA_early <- ggplot(df_min_early, aes(x = DAYS_TX_TO_RESULT, y =factor(PAT_ID), color = Status)) +
  geom_point(size = 3) +
  labs( title = "Days till first DSA result (Early)") +
  theme_minimal() +   
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60 days

###################################################
# Ridgeplots: Serotype distribution across days 
###################################################

# Create a ridge plot for all early patients
ggplot(df_unique_early, aes(x = DAYS_TX_TO_RESULT, y = type, height = ..density..)) +
  geom_density_ridges(stat = "density", scale = 2, fill = "lightblue", color = "blue") +
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution Across Days (Early patient group)") +
  theme_minimal()  +   
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1)  # Add vertical line at 60 days

# get patient status column
df_unique_early2 <- early_df_separated %>%
  distinct(DAYS_TX_TO_RESULT, SEROTYPE, Status)
df_unique_early2$type<- extract_letters(df_unique_early2$SEROTYPE)
df_unique_early2$Status <- factor(df_unique_early2$Status)


# NR/R patients overlap
early_ridge <- ggplot(df_unique_early2, aes(x = DAYS_TX_TO_RESULT,fill=Status, y = type, height = ..density..)) + 
  geom_density_ridges(stat = "density", scale = 2, alpha = 0.5) + 
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution across time \n (Early group)") + 
  theme_minimal() + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1) 


# NR/R patients overlap
ggplot(df_unique_early2, aes(x = DAYS_TX_TO_RESULT, y = type, height = ..density..)) + 
  geom_density_ridges(stat = "density", scale = 2) + 
  labs(x = "Days to Result", y = "Serotype", title = "Serotype Distribution Across Days (Early group)") + 
  theme_minimal() + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~Status)

plot_grid(firstDSA_early, early_ridge, time_early, time_DSA_early, time_ABMR_early)

#######
# cox analysis
#####
# if a serotype is stored as a list, will expand this list into additional rows
hist_separated <- hist %>%
  separate_rows(SEROTYPE, sep = ",") %>%
  drop_na(SEROTYPE)

hist_separated$event <- ifelse(is.na(hist_separated$SEROTYPE),0,1)
# Install and load the survival package if not already installed
install.packages("survival")  # Uncomment this line if you need to install the package
library(survival)

# Example dataset structure
# Ensure your dataset contains the following columns:
#  - time: Survival time (e.g., days, months, etc.)
#  - status: Censoring status (1 = event occurred, 0 = censored)
#  - covariates: Predictor variables (e.g., age, treatment group, etc.)

# Fit Cox proportional hazards model
cox_model <- coxph(Surv(DAYS_TX_TO_RESULT, event) ~ Status,DSA_CLASS, data = hist_separated)

# Summarize the Cox model
summary(cox_model)
