#!/bin/sh
#SBATCH --job-name=nWF_popSize_dynasty
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim

Rscript  burnIn_False.R
Rscript  checkFileNWF2.R

slim dynasty_nWF_populationsize.slim
wait
