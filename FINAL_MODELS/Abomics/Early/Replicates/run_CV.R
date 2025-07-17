source("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/SLIDE_CV_EARLY.R")
library(SLIDE)
library(EssReg)
yaml_path = "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/slide_cv_2.yaml"
SLIDEcv_f(yaml_path, nrep = 100, k = 18)