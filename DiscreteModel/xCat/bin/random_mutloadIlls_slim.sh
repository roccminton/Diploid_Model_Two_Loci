#!/bin/sh
#SBATCH --job-name=mutload_random
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim

Rscript  checkFile.R

slim random_mutloadIlls.slim
wait

