#!/bin/bash
# T9 (b): shell psi oscillation factor study. usage: sbatch --export=ALL,Z=80,NCH=8,CD=20 scripts/slurm/l5_t9_osc.sbatch
#SBATCH --ntasks=1
#SBATCH --partition=sequana_cpu
#SBATCH --time=02:00:00
#SBATCH --job-name=t9-osc
#SBATCH --output=/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/l5_oneway/t9d/slurm-osc-%j.out
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/t9osc-$SLURM_JOB_ID; mkdir -p $TMPDIR
scripts/aeroenv.sh env PYTHONPATH=$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python tools/l5_t9_psi_osc.py run --z ${Z} --nchord ${NCH} --cd ${CD} --out logs/check/l5_oneway/t9d/osc/z${Z}_n${NCH}_cd${CD}
