#!/bin/bash
# usage: sbatch -J <tag> -p <partition> --time=<t> job.sbatch <tree T7e|T14m> <case> <es>
# PREDICTION (written before submission; AE = our model, gap = AE/CCX-1, W gap per T21 = -9.6 % T7e, -20.5 % T14m for tip_axial_tension):
#  - es 0.25 tip_axial_tension: reproduces T21 W gaps (T7e ~ -9.6 %, T14m ~ -20.5 %).
#  - tip_flap / tip_edge / tip_torsion: if the cause is the eccentric-load amplifier, each of these W gaps is < 5 % on both trees
#    (flap/edge as in the docstring, 1.7 / 0.5 % on the older tree); torsion is the least certain (GJ moved 5 -> 10 %), expect 5-12 %.
#    If instead flap or edge or torsion shows ~ -10..-20 % on T14m, the formulation/mesh (warped TE strips) is the driver.
#  - es 0.125 tip_axial_tension T14m: gap falls toward <= 10 % (docstring trend 12.6 -> 3.3 % with refinement); if it stays ~ -20 % the gap is a model difference.
#  Runtimes: es 0.25 ~ 20-60 min (T21 timed out the 4-case run at 90 min, one case ran in under 90); es 0.125 ~ 3-6 h.
#SBATCH --ntasks=8
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t31/pc-%x-%j.log
source /scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh
W=$SCRATCH/wt/input-fix-$1
export TMPDIR=$SCRATCH/tmp/t31-$SLURM_JOB_ID; mkdir -p $TMPDIR
export PYTHONPATH=$W.pp:$W
cd $W
python /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/t31/percase.py $W $2 $3
