#!/bin/bash
# usage: sbatch ... reg.sbatch <group> <scope-file>   (regression without --write, JSON)
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/sp1-$SLURM_JOB_ID; mkdir -p $TMPDIR
echo "host $(hostname) rev $(git rev-parse --short HEAD) start $(date +%T)"
/usr/bin/time -f "elapsed %e s maxrss %M kB" scripts/aeroenv.sh python -B tools/validation_matrix.py regression --group $1 --scope "$(cat $2)" --json > logs/check/l3/store-pass1/reg-$1.json 2> logs/check/l3/store-pass1/reg-$1.err
echo "exit $? $(date +%T)"; tail -n 2 logs/check/l3/store-pass1/reg-$1.err
