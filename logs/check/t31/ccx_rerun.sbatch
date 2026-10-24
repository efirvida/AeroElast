#!/bin/bash
# Re-run CCX on the preserved es 0.125 deck of 11623429 (rc 16, stdout lost) to read its diagnostic.
#SBATCH --ntasks=8
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t31/ccx-rerun-%j.log
source /scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh
export OMP_NUM_THREADS=8 CCX_NPROC_EQUATION_SOLVER=8
cd $SCRATCH/tmp/t31-11623429/tmpl4q4d8aj
/usr/bin/time -v ccx b; echo "ccx exit $?"
