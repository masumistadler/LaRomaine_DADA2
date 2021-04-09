#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=0-30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --job-name="track_reads"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load gcc/9.3.0 r/4.0.2 openmpi/4.0.3

#command to execute
Rscript ~/projects/def-pauldel/Jobs/TrackReads_slurm.R
#---------------------------------------------
