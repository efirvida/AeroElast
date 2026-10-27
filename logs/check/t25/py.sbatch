#!/bin/bash
# usage: sbatch -J <tag> -p <part> --time=<t> py.sbatch <script.py>   (runs in the pinned a35bba3 worktree)
#SBATCH --ntasks=4
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t25/%x-%j.log
W=${W:-/scratch/leahk/eduardo.donestevez/wt/full-a35bba3}
export PYTHONPATH=$W.pp:$W PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
export TMPDIR=$SCRATCH/tmp/t25-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd $W || exit 2
/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh python -W ignore "$1" 2>&1 | grep -vE "^\s+(Step|Generat|Creating|Dedup|Load|Assign|Built|Blade mesh)"
