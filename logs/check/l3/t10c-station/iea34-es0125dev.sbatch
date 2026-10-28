#!/bin/bash
#SBATCH --job-name=t10cs-iea34-0125dev
#SBATCH --partition=sequana_cpu_dev
#SBATCH --time=00:20:00
#SBATCH --ntasks=4
#SBATCH --output=/scratch/leahk/eduardo.donestevez/wt/l3-pass2/logs/check/l3/t10c-station/iea34-es0125dev-%j.log
set -euo pipefail
cd /scratch/leahk/eduardo.donestevez/wt/l3-pass2
export TMPDIR=$SCRATCH/tmp/t10cs-$SLURM_JOB_ID; mkdir -p $TMPDIR
export PYTHONPATH=/scratch/leahk/eduardo.donestevez/wt/l3-pass2/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
echo "worktree rev $(cat logs/check/l3/t10c-station/revision)"
/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh python -u tools/s2_ring_metric_blade.py --yaml tests/IEA-3.4-130-RWT.yaml --element-size 0.125 --out logs/check/l3/t10c-station/iea34-es0125dev.json
