library(readxl)
library(xlsx)
library(dplyr)
library(stringr)

# helper functions
remove_underscore<- function(string){
  new_string <- sub("_.*$", "", string)
  return(new_string)
}

remove_ending_number_pattern <- function(s) {
  # Remove the pattern ".[numbers]." at the end of the string
  sub("\\.[0-9]+\\.$", "", s)
}

remove_prefix <- function(string) {
  return(sub("^[^_]*_", "", string))
}

################################
# cleaning original spreadsheet
################################

# antigens <- read_excel("ABomics/pitt_cohort/antigens/One Lambda HLA Sequences.xlsx", sheet=2)
# 
# 
# # fixing syntax issues 
# ant <- sapply(colnames(antigens) , function(s){ gsub(" ", "", s)})
# ant2 <- sapply(ant, function(s){ gsub("-", ".", s)})
# ant3 <- sapply(ant2 , function(s){ gsub("ClassMICA", "MICA", s)})
# colnames(antigens)<- ant3


###########################################################################
# checking antigen types and classes 
###########################################################################

antigens_serotype_mapping_reference <- readRDS("/ix/djishnu/Marisa/ABomics/pitt_cohort/antigens/antigens_serotype_mapping_reference.rds")
all_antigens <- na.omit(unique(unlist(antigens_serotype_mapping_reference)))


check_RA <- function(string) {
  if (startsWith(string, "A") | startsWith(string, "B") | startsWith(string, "C")) {
    return("Antigen")
  } 
  else if (startsWith(string, "D")) {
    return("Receptor")
  } else {
    return("NonDSA")
  }
}


classtype<- function(string) {
  if (startsWith(string, "A") | startsWith(string, "B") | startsWith(string, "C")) {
    return("I")
  } else if (startsWith(string, "D")) {
    return("II")
  } else {
    return("MICA")
  }
}

aa <- data.frame("Name" = all_antigens)
names(all_antigens) <- all_antigens
aa$type <- sapply(all_antigens, check_RA)
aa$class <- sapply(all_antigens, classtype)


