#!/bin/sh
#SBATCH --job-name=nWF_popsize_random
#SBATCH --partition=batch
#SBATCH --time=12:00:00
module load R/3.5.1
module load slim
Rscript  checkLambda.R

slim random_popsize_smallGenome_difLambda.slim
wait
