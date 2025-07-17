#!/bin/bash

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH --job-name=ER
#SBATCH --output=/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/output.txt
#SBATCH --error=/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/error.txt
#SBATCH --mail-user=mia100@pitt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=100GB
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16

module purge
module load gcc/12.2.0 r/4.3.0
rscript="/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/run_CV.R"
rout= "/ix/djishnu/Marisa/ABomics/FINAL_MODELS/Abomics/Early/Replicates/run_CV.Rout"

# Run the R script
R CMD BATCH $rscript $rout

