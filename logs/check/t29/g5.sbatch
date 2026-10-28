#!/bin/bash
# T29 row 5: the mesh-convergence test with per-CCX-mode gaps, printed (-s), on the main repo HEAD.
# PREDICTION (before the run): every CCX mode whose 0.5 m gap is > 0.2 % falls 0.5 -> 0.35 m, so the
# convergence loop passes and both residual lines print (row 5 maps). Low worst < 1 %. High worst at
# 0.25 m about 5.2 % (docstring 5.23 %), so the high-mode residual assertion is expected to FAIL against
# 5 %: a recorded discretisation finding, not a bound to move. If the loop still fails on a CCX mode,
# that mode really does not converge 0.5 -> 0.35 (finding).
#SBATCH --ntasks=8
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t29/g5-%j.log
R=/scratch/leahk/eduardo.donestevez/fem-shell
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR=$SCRATCH/tmp/t29-$SLURM_JOB_ID; mkdir -p $TMPDIR
cd $R && echo "rev=$(cat .git/HEAD 2>/dev/null)"
scripts/aeroenv.sh python -m pytest tests/validation/parity/test_blade_iea15mw_mesh_convergence.py -s -q -rf --tb=short -p no:cacheprovider -W ignore --basetemp=$TMPDIR/bt -o addopts="" 2>&1 | grep -vE "Warning|warnings.warn"
rm -rf $TMPDIR
