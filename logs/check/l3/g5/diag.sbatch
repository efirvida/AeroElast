#!/bin/bash
#SBATCH --job-name=g5-mac
#SBATCH --partition=sequana_cpu
#SBATCH --time=01:30:00
#SBATCH --ntasks=8
#SBATCH --output=/scratch/leahk/eduardo.donestevez/wt/l3-g5/logs/check/l3/g5/diag-%j.log
set -euo pipefail
W=/scratch/leahk/eduardo.donestevez/wt/l3-g5; cd $W
export TMPDIR=$SCRATCH/tmp/g5-mac-$SLURM_JOB_ID; mkdir -p $TMPDIR
export PYTHONPATH=$W/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "worktree rev $(cat logs/check/l3/g5/revision) (+ uncommitted T12 diff)"
/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh python -u tools/diagnose_g5_mode7.py --meshes 0.5 0.35 0.25 --modes 10 --out logs/check/l3/g5/mac.json 2>&1 | grep -v -i "warning"
rm -rf $TMPDIR
