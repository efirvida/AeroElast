#!/bin/bash
#SBATCH --partition=sequana_cpu
#SBATCH --time=90
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/check/t10/%x-%j.log
set -eu
cd /scratch/leahk/eduardo.donestevez/fem-shell
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
scripts/aeroenv.sh python -B tools/t10_ccx_arbiter.py "$1" --out "logs/check/t10/$1"
