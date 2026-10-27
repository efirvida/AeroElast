#!/bin/bash
#SBATCH --job-name=t29-g5
#SBATCH --partition=sequana_cpu
#SBATCH --ntasks=4
#SBATCH --time=02:00:00
#SBATCH --output=logs/check/t29/g5-%j.out
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/t29-$SLURM_JOB_ID; mkdir -p $TMPDIR
export RAYON_NUM_THREADS=4
scripts/aeroenv.sh python tools/validation_matrix.py regression --group 5 --write --json > logs/check/t29/g5.json 2> logs/check/t29/g5.err
echo exit $?
