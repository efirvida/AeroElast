#!/bin/bash
#SBATCH -p sequana_cpu_dev
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH -J l5_twoway_of
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/l5_twoway/%x-%j.out
# NOT SUBMITTED. flex ~0.11 sim-s/CPU-s (60 s ~ 9 min), rigid ~0.36 (3 min). Serial, so the smallest request:
#   sbatch -p sequana_cpu_dev --ntasks=1 --time=00:20:00 logs/check/l5_twoway/run-openfast.sh
set -uo pipefail
OF=/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast
DECKS="${SCRATCH:-/scratch/leahk/eduardo.donestevez}/l5_twoway/openfast"
for f in flex rigid; do
  ( cd "$DECKS/$f" && "$OF" run.fst > run.log 2>&1; echo "$f rc=$? $(tail -n 3 run.log | tr '\n' ' ')" )
done
