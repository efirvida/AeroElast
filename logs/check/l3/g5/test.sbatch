#!/bin/bash
#SBATCH --job-name=g5-test
#SBATCH --partition=sequana_cpu_bigmem
#SBATCH --time=01:30:00
#SBATCH --ntasks=8
#SBATCH --output=/scratch/leahk/eduardo.donestevez/wt/l3-g5/logs/check/l3/g5/test-%j.log
W=/scratch/leahk/eduardo.donestevez/wt/l3-g5; cd $W
export TMPDIR=$SCRATCH/tmp/g5-test-$SLURM_JOB_ID; mkdir -p $TMPDIR
export PYTHONPATH=$W/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "worktree rev $(cat logs/check/l3/g5/revision) (+ uncommitted T12 diff)"
/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh python -m pytest tests/validation/parity/test_blade_iea15mw_mesh_convergence.py -s -q -rf --tb=short -p no:cacheprovider -W ignore --basetemp=$TMPDIR/bt -o addopts="" 2>&1 | grep -vE "Warning|warnings.warn"
rm -rf $TMPDIR
