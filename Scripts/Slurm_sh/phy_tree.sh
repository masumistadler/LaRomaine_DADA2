#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=3-00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=24
#SBATCH --job-name="phy.tree"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#creating a variable for the blast database

#to load the executables
module load gcc/9.3.0 r/4.2.2 openmpi/4.0.3

#command to execute
Rscript /home/mstadler/projects/def-pauldel/mstadler/Jobs/8_phy.tree_slurm.R
#choose --save if you want #to #save your results on the R interface
