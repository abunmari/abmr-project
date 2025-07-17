setwd("~/JishnuLab/ABomics_cleaned")


source("~/JishnuLab/ABomics_cleaned/Analyses/cross_preds/cp_SLIDE.R")

# TESTING CROSS PREDICTION PERFORMANCE AFTER REMOVE of LFs
# LATE BASED ON COHORT SPECIFIC MEDIAN
late_cp_yaml = "FINAL_MODELS/Abomics/Late/input.yaml"
canada_late_x = "canada_cohort/Data/AB/late_x.csv"
canada_late_y =  "canada_cohort/Data/AB/late_y.csv"

#pdf("~/JishnuLab/ABomics_7_10/FINAL_MODELS/Abomics/Late/FIGURES/cross_pred_cohort_median.pdf")

ab_late_marg <- load_data_cp(yaml_path = late_cp_yaml, 
                             val_x = canada_late_x,
                             val_y =canada_late_y,
                             interactions = FALSE, scale = TRUE)



l1 <- load_data_cp(yaml_path = late_cp_yaml,  # train auc = 0.93, val_auc = 0.83
                  val_x = canada_late_x,
                  val_y =canada_late_y,
                  interactions = F, scale = TRUE, marg_z = c(1,2,7,10))

# remove Z10
l2 <- load_data_cp(yaml_path = late_cp_yaml, # train_auc = 0.9, val_auc = 0.9
                  val_x = canada_late_x,
                  val_y =canada_late_y,
                  interactions = F, scale = TRUE, marg_z = c(1,2,7))
# Remove Z7
l3<- load_data_cp(yaml_path = late_cp_yaml, # train_auc = 0.81, val_auc = 0.73
                  val_x = canada_late_x,
                  val_y =canada_late_y,
                  interactions = F, scale = TRUE, marg_z = c(1,2,10))
# Remove Z2
l4<- load_data_cp(yaml_path = late_cp_yaml, # train_auc = 0.89, val_auc = 0.84
                  val_x = canada_late_x,
                  val_y =canada_late_y,
                  interactions = F, scale = TRUE, marg_z = c(1,7,10))

# Remove Z1
l5<- load_data_cp(yaml_path = late_cp_yaml, # train_auc = 0.85, val_auc = 0.97
                  val_x = canada_late_x,
                  val_y =canada_late_y,
                  interactions = F, scale = TRUE, marg_z = c(2,7,10))
###############################################################################
# early model on late model ---------------------------------------------------
###############################################################################

# setwd("~/JishnuLab/ABomics_cleaned")
early_cp_yaml = "FINAL_MODELS/Abomics/Early/input.yaml"


pitt_early_x = "FINAL_MODELS/Abomics/Early/x.csv"
pitt_early_y = "FINAL_MODELS/Abomics/Early/y.csv"

pitt_late_x = "FINAL_MODELS/Abomics/Late/x.csv"
pitt_late_y = "FINAL_MODELS/Abomics/Late/y.csv"

#canada_early_x = "canada_cohort/Data/AB/early_x.csv"
#canada_early_y =  "canada_cohort/Data/AB/early_y.csv"

# early signature on late signature - works
early_late <- load_data_cp(yaml_path = early_cp_yaml, 
                        val_x = pitt_late_x,
                        val_y = pitt_late_y,
                        interactions = FALSE, scale = TRUE)

# late signature on early - doesn't work
late_early <- load_data_cp(yaml_path = late_cp_yaml, 
                           val_x = pitt_early_x,
                           val_y = pitt_early_y,
                           interactions = FALSE, scale = TRUE)

# late signature on early - doesn't work
late_early <- load_data_cp(yaml_path = late_cp_yaml, 
                           val_x = pitt_early_x,
                           val_y = pitt_early_y,
                           interactions = FALSE, scale = TRUE, marg_z = c(2,7))


# testing if there is a pair of late lfs that can predict early 
pairwise_combinations <- function(vec) {# Use combn to generate pairwise combinations
  combinations <- combn(vec, 2, simplify = FALSE)
  return(combinations)
}

late_lfs <- c(1,2,7,10)
lz = pairwise_combinations(late_lfs)
perf_list <- NULL
for(i in 1:length(lz)){
  zpair = lz[[i]]
late_early <- load_data_cp(yaml_path = late_cp_yaml, 
                             val_x = pitt_early_x,
                             val_y = pitt_early_y,
                             interactions = FALSE, scale = TRUE, marg_z = zpair)
                             perf_list[[i]] <- paste0(zpair[1],"_",zpair[2],  " train_auc: " ,late_early$train_auc," val_auc: " , late_early$val_auc)
}
# LF 2 and 7 from late model can predict on early patients
late_early_Z_2_7<- load_data_cp(yaml_path = late_cp_yaml, 
                           val_x = pitt_early_x,
                           val_y = pitt_early_y,
                           interactions = FALSE, scale = TRUE, marg_z = c(2,7))

############################################################################
# now testing signatures on toronto cohort
############################################################################

pitt_early_x = "FINAL_MODELS/Abomics/Early/x.csv"
pitt_early_y = "FINAL_MODELS/Abomics/Early/y.csv"
canada_late_x = "canada_cohort/Data/AB/late_x.csv"
canada_late_y =  "canada_cohort/Data/AB/late_y.csv"
canada_early_x = "canada_cohort/Data/AB/early_x.csv"
canada_early_y =  "canada_cohort/Data/AB/early_y.csv"

early_late_tor <- load_data_cp(yaml_path = early_cp_yaml, 
                           val_x =canada_late_x,
                           val_y =canada_late_y,
                           interactions = FALSE, scale = TRUE)

late_early_tor  <- load_data_cp(yaml_path = late_cp_yaml, 
                                val_x =canada_early_x,
                                val_y =canada_early_y,
                                interactions = FALSE, scale = TRUE, labels= c("Pitt Late", "Toronto early"))
############################################################################
# now testing signatures on combined late cohort (pitt + toronto )
############################################################################
late_cp_yaml = "FINAL_MODELS/Abomics/Late/input.yaml"

late_cp_yaml = "FINAL_MODELS/Abomics/Late/input.yaml"

all_late_x <-"Analyses/cross_preds/all_late_x.csv"
all_late_y <-"Analyses/cross_preds/all_late_y.csv"

late_tor_pitt<- load_data_cp(yaml_path = late_cp_yaml, 
                               val_x =all_late_x,
                               val_y =all_late_y,
                               interactions = FALSE, scale = TRUE,
                             labels = c("Pitt Late", "Pitt & Toronto Late"))

