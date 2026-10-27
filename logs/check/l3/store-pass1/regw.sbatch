#!/bin/bash
# usage: sbatch ... regw.sbatch <group>   (regression --write over the group's own source files)
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/sp1-$SLURM_JOB_ID; mkdir -p $TMPDIR
echo "host $(hostname) rev $(git rev-parse --short HEAD 2>/dev/null) start $(date +%T)"
/usr/bin/time -f "elapsed %e s maxrss %M kB" scripts/aeroenv.sh python -B tools/validation_matrix.py regression --group $1 --write --json > logs/check/l3/store-pass1/regw-$1.json 2> logs/check/l3/store-pass1/regw-$1.err
echo "exit $? $(date +%T)"; tail -n 3 logs/check/l3/store-pass1/regw-$1.err
