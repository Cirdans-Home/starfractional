#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --job-name=slanczos
#SBATCH -N 1
#SBATCH --cpus-per-task=256
#SBATCH --time 10-00:00:00

export OMP_NUM_THREADS=256
module load matlab/R2021a
matlab -batch "autonomous_solution_table"

