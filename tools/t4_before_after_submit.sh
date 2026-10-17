#!/bin/bash
# Submit the T4 before/after measurement matrix (oml vs mid) as SLURM jobs.
# Usage: tools/t4_before_after_submit.sh [--go]   (default: print the sbatch commands only)
# Same command on both sides: only BladeMesh's default skin_reference is forced
# (tools/with_skin_reference.py). Logs: logs/check/t4/<target>-<ref>-<jobid>.log
cd "$(dirname "$0")/.." || exit 1
mkdir -p logs/check/t4
GO=0; [ "$1" = "--go" ] && GO=1
run() {  # name ref minutes ntasks command...
  local name=$1 ref=$2 minutes=$3 ntasks=$4; shift 4
  local script=logs/check/t4/$name-$ref.sbatch
  cat > "$script" <<EOF
#!/bin/bash
#SBATCH --partition=sequana_cpu
#SBATCH --time=$minutes
#SBATCH --ntasks=$ntasks
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/check/t4/$name-$ref-%j.log
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RAYON_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
scripts/aeroenv.sh python -B tools/with_skin_reference.py $ref $*
EOF
  if [ $GO = 1 ]; then sbatch "$script"; else echo "sbatch $script   # $*"; fi
}
for ref in oml mid; do
  run s1    $ref 40  1 tools/s1_reference_sweep.py $ref 0.25
  run s2    $ref 120 1 -m pytest -o addopts= -q -rA -p no:cacheprovider tests/test_iea15mw_s2_static_prescribed.py
  run s4    $ref 60  1 -m pytest -o addopts= -q -rA -p no:cacheprovider tests/test_iea15mw_s4_rotating_modal.py
  run s4c   $ref 180 12 tools/s4_modal_convergence.py 0.25 0.125 --rpm=7.56
  run v02   $ref 60  1 -m pytest -o addopts= -q -rA -p no:cacheprovider tests/validation/blade/test_iea15mw_v02_natural_frequencies.py
  run v03   $ref 60  1 -m pytest -o addopts= -q -rA -p no:cacheprovider tests/test_iea15mw_v03_static_gravity.py
done
