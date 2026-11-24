#!/bin/bash
#SBATCH -J l5-s6
#SBATCH -p sequana_cpu_bigmem
#SBATCH --ntasks=4
#SBATCH --time=60:00
#SBATCH -o logs/check/l5_oneway/l5-s6-%j.log
export TMPDIR=$SCRATCH/tmp/l5-s6-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd /scratch/leahk/eduardo.donestevez/fem-shell
scripts/aeroenv.sh python -m pytest -q -p no:cacheprovider --basetemp=$TMPDIR/pt tests/test_iea15mw_s6_oneway_dynamic.py tests/test_iea15mw_s5_oneway.py -m slow -W ignore 2>&1 | tail -30
