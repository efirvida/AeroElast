#!/bin/bash
# usage: sbatch -J <tag> job <tree T7e|T14m> <frd> <es>
# PREDICTION (before run): if dtheta/dz AE/CCX-1 is ~ -5 % or less in mid-span bands while the tip theta gap is ~ -23..-26 %,
# the torsion W gap is loaded-end distortion (ill-posed load), not GJ; if mid-span bands also show ~ -20 %, it is a GJ gap.
#SBATCH --ntasks=8
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t31/station-%x-%j.log
source /scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh
W=$SCRATCH/wt/input-fix-$1
export TMPDIR=$SCRATCH/tmp/t31-$SLURM_JOB_ID; mkdir -p $TMPDIR
export PYTHONPATH=$W.pp:$W
cd $W
python /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t31/station_twist.py $W $2 $3 $4
