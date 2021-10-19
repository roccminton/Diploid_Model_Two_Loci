#!/bin/sh
#SBATCH --job-name=decline_dynasty
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1
module load slim

Rscript  checkFile2.R

slim dynasty_decline.slim
wait
