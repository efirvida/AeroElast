#!/bin/bash
#SBATCH --job-name=l3-g23
#SBATCH --ntasks=8
#SBATCH --time=01:30:00
#SBATCH --partition=sequana_cpu_bigmem
#SBATCH --output=/scratch/leahk/eduardo.donestevez/wt/l3-g23/logs/check/l3/g23/job-%j.out
set -u
export TMPDIR=$SCRATCH/tmp/g23-$SLURM_JOB_ID
mkdir -p "$TMPDIR"
cd /scratch/leahk/eduardo.donestevez/wt/l3-g23
export PYTHONPATH=$PWD/src
/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh python -m pytest \
    tests/validation/blade/test_blade_iea15mw_validation.py -s -v -p no:cacheprovider \
    --basetemp="$TMPDIR/pytest" 2>&1 | tail -n 400
