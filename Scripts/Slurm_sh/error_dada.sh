#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=10-00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=24 #34
#SBATCH --job-name="error_dada"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load gcc/9.3.0 r/4.0.2 openmpi/4.0.3

#command to execute
Rscript ~/projects/def-pauldel/mstadler/Jobs/3_learnError_dada_slurm.R
#---------------------------------------------
