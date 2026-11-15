#!/bin/bash
#SBATCH -J l5-ad
#SBATCH -p sequana_cpu_dev
#SBATCH --ntasks=4
#SBATCH --time=19:00
#SBATCH -o logs/check/l5_oneway/l5-ad-%j.log
export TMPDIR=$SCRATCH/tmp/l5-ad-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd /scratch/leahk/eduardo.donestevez/fem-shell
scripts/aeroenv.sh python tools/l5_oneway_bd_parity.py --loads ad --out logs/check/l5_oneway/ad.json --work logs/check/l5_oneway/bd_ad
