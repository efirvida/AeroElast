#!/bin/bash
#SBATCH -J l5-bem
#SBATCH -p sequana_cpu
#SBATCH --ntasks=4
#SBATCH --time=40:00
#SBATCH -o logs/check/l5_oneway/l5-bem-%j.log
export TMPDIR=$SCRATCH/tmp/l5-bem-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd /scratch/leahk/eduardo.donestevez/fem-shell
scripts/aeroenv.sh python tools/l5_oneway_bd_parity.py --loads bem --out logs/check/l5_oneway/bem.json --work logs/check/l5_oneway/bd_bem
