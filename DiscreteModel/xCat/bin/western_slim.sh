#!/bin/sh
#SBATCH --job-name=western_slim
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim
Rscript  checkFile.R

slim western_xcat.slim
wait
