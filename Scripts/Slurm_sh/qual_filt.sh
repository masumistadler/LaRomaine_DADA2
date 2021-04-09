#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=1-00:00
#SBATCH --mem=10G #10G
#SBATCH --cpus-per-task=12 #48, max 48
#SBATCH --job-name="qual_filt_dada"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load gcc/9.3.0 r/4.0.2 openmpi/4.0.3

#command to execute
Rscript ~/projects/def-pauldel/mstadler/Jobs/2_quality_filt_slurm.R
#---------------------------------------------
