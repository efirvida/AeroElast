#!/bin/bash
# T9 (a): z45 solid torsion length convergence. usage: sbatch --export=ALL,NCH=8 scripts/slurm/l5_t9_solid.sbatch
# 11744952-55 (ntasks 3, 24 GB) were OOM-killed: three ccx load cases run at once; ntasks 12 gives 96 GB.
#SBATCH --ntasks=12
#SBATCH --partition=sequana_cpu
#SBATCH --time=04:00:00
#SBATCH --job-name=t9-solid
#SBATCH --output=/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/l5_oneway/t9d/solid/slurm-%j.out
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/t9solid-$SLURM_JOB_ID; mkdir -p $TMPDIR
OUT=logs/check/l5_oneway/t9d/solid/z45_A_n${NCH}
mkdir -p $OUT
scripts/aeroenv.sh env PYTHONPATH=$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python tools/l5_t9_arbiter.py solid --z 45 --level A --nchord ${NCH} --out $OUT
