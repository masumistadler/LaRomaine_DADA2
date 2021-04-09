#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=0-60:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=3
#SBATCH --job-name="Install_Rpackages"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load StdEnv/2020 gcc/9.3.0 r/4.0.2 openmpi/4.0.3

#command to execute
Rscript /home/mstadler/projects/def-pauldel/mstadler/Jobs/install_packages.R
#---------------------------------------------
