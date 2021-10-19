#!/bin/sh
#SBATCH --job-name=nWF_popsize_random
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim
Rscript  burnIn_True.R
Rscript  checkFileNWF.R

slim random_nWF_populationsize.slim
wait
