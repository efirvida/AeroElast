#!/bin/bash
#SBATCH --job-name=restamp-d5d0fd2
#SBATCH --output=/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/restamp-d5d0fd2/batch-%x-%j.log
# Restamp of every blade-dependent store group at the blade input audit T9 tree (a35bba3 = d5d0fd2 code),
# odd/tasks/blade-input-audit.md Handoff 4 item 1. Pinned worktree + shim; --write lands in the worktree's
# docs/validation and is reviewed by git diff before being carried to the branch.
# Prediction stated before running: blade groups move (geometry changed); G11 (BEM, aero only) same.
S=/scratch/leahk/eduardo.donestevez; R=$S/fem-shell; W=$S/wt/restamp-a35bba3; B=$R/logs/check/restamp-d5d0fd2
cd $W || exit 2
export PATH=$S/gitshim-d5d0fd2:$PATH PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH=$S/wt/restamp-a35bba3.pp
scripts/aeroenv.sh python -c "import aeroelast,sys; sys.exit(0 if aeroelast.__file__.startswith('$S/wt/restamp-a35bba3.pp/') else 3)" || { echo "aeroelast not from the shim"; exit 3; }
echo "job=$SLURM_JOB_ID host=$HOSTNAME wt=$(/usr/bin/git rev-parse --short HEAD) stamp=$(git rev-parse --short HEAD)"
reg() { tag=$1; shift
  export TMPDIR=$S/tmp/restamp-$SLURM_JOB_ID-$tag RAYON_NUM_THREADS=2; mkdir -p $TMPDIR
  for g in "$@"; do t0=$(date +%s); sc=""; [ "$g" = 43 ] && sc="--scope tests/test_iea15mw_s3_newmark_transient.py"
    scripts/aeroenv.sh python tools/validation_matrix.py regression --group $g $sc --write --json > $B/g$g.json 2> $B/g$g.err
    echo "$tag group $g exit $? $(( $(date +%s) - t0 )) s"; done; }
cap() { tag=$1; shift
  export TMPDIR=$S/tmp/restamp-$SLURM_JOB_ID-cap$tag RAYON_NUM_THREADS=4 PYTHONPATH=$S/wt/restamp-a35bba3.pp:$R/logs/check/b3-1r B3_CAPTURE=$B/capture-$tag.pkl; mkdir -p $TMPDIR
  t0=$(date +%s)
  scripts/aeroenv.sh python -m pytest -q -s -p no:cacheprovider -p capture_plugin --basetemp=$TMPDIR/pt "$@" -rA > $B/pytest-$tag.log 2>&1
  echo "capture $tag exit $? $(( $(date +%s) - t0 )) s"; tail -n 1 $B/pytest-$tag.log; }
reg s1 58 11 > $B/s1.log 2>&1 &
reg s2 20 32 34 > $B/s2.log 2>&1 &
reg s3 5 23 > $B/s3.log 2>&1 &
reg s4 27 43 > $B/s4.log 2>&1 &
cap s4 tests/test_iea15mw_s4_rotating_modal.py > $B/cap-s4.log 2>&1 &
cap v05 tests/test_iea15mw_v05_structural_properties.py tests/test_iea15mw_s6_oneway_dynamic.py > $B/cap-v05.log 2>&1 &
wait; cat $B/s?.log $B/cap-*.log; cd $W && /usr/bin/git diff --stat -- docs/validation
