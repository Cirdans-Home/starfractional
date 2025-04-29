#!/bin/bash
#SBATCH --partition=cl2
#SBATCH --job-name=slanczos
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --time 10-00:00:00

export OMP_NUM_THREADS=48
module load matlab/R2021a
matlab -batch "autonomous_solution_table"

