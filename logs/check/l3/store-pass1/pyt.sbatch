#!/bin/bash
# usage: sbatch -J name -p P --ntasks=4 --time=T -o logs/check/l3/store-pass1/%x-%j.log logs/check/l3/store-pass1/pyt.sbatch <out-tag> <pytest args...>
cd /scratch/leahk/eduardo.donestevez/fem-shell
export TMPDIR=$SCRATCH/tmp/sp1-$SLURM_JOB_ID; mkdir -p $TMPDIR
tag=$1; shift
echo "host $(hostname) rev $(git rev-parse --short HEAD) start $(date +%T)"
PYTHONPATH=tests /usr/bin/time -f "elapsed %e s maxrss %M kB" scripts/aeroenv.sh python -B -m pytest -o addopts= -s -v -rA -p no:cacheprovider --basetemp=$TMPDIR/bt "$@" > logs/check/l3/store-pass1/pyt-$tag.txt 2>&1
echo "exit $? $(date +%T)"; tail -n 3 logs/check/l3/store-pass1/pyt-$tag.txt
