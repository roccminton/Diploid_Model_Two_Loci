#!/bin/sh
#SBATCH --job-name=createGenes
#SBATCH --partition=batch
#SBATCH --time=10:00:00
module load R/3.5.1

Rscript  geneErstellen2.R
wait
