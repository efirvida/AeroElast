#!/bin/bash
# usage: sbatch -J <tag> -p <part> --time=<t> rerun.sbatch <idfile>
# T25: re-run the failures of full suite 11623361 in the same pinned worktree (a35bba3) to read their messages.
#SBATCH --ntasks=8
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t25/rerun-%x-%j.log
W=/scratch/leahk/eduardo.donestevez/wt/full-a35bba3
export PYTHONPATH=$W.pp PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR=$SCRATCH/tmp/t25-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd $W || exit 2
echo "job=$SLURM_JOB_ID rev=$(git rev-parse --short HEAD)"
scripts/aeroenv.sh python -m pytest $(cat "$1") -q -rf --tb=short -p no:cacheprovider -W ignore --basetemp=$TMPDIR/bt -o addopts="" 2>&1 | grep -vE "Warning|warnings.warn"
