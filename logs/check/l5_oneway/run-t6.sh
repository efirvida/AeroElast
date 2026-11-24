#!/bin/bash
#SBATCH -p sequana_cpu
#SBATCH --ntasks=8
#SBATCH -J l5_oneway_rated
#SBATCH -o /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/l5_oneway/%x-%j.out
#SBATCH -e /scratch/leahk/eduardo.donestevez/fem-shell/logs/check/l5_oneway/%x-%j.err
#SBATCH --time=04:00:00
#
# L5 T6 (odd/tasks/l5-one-way.md): S-5 rated run with the production code path, one-way.
#   BEM (aeroelast-bem-fsi) + production rotor shell solver (aeroelast) over preCICE, no deformation feedback
#   (deformed_twist / deformed_radius false), tilt -6 / precone -4, 10.59 m/s, 7.56 rpm, gmsh N20 AR 1 mid mesh.
# Usage:  sbatch logs/check/l5_oneway/run-t6.sh            # full 100 s
#         sbatch -p sequana_cpu_dev --time=00:20:00 logs/check/l5_oneway/run-t6.sh smoke   # 0.5 s timing smoke
# Output: $SCRATCH/l5_oneway/<mode>-<jobid>/  (results/, bem_fsi_results/, precice-run/, *.log, watchpoint files)
set -uo pipefail

MODE="${1:-full}"
REPO=/scratch/leahk/eduardo.donestevez/fem-shell
CASE="$REPO/cases/l5_oneway_rated"
case "$MODE" in
  full)  FLUID_YAML="$CASE/fluid/fluid.yaml"; SOLID_YAML="$CASE/solid/solid.yaml" ;;
  smoke) FLUID_YAML="$CASE/smoke/fluid.yaml"; SOLID_YAML="$CASE/smoke/solid.yaml" ;;
  *) echo "mode must be full|smoke"; exit 2 ;;
esac

RUN="$SCRATCH/l5_oneway/${MODE}-${SLURM_JOB_ID:-local}"
mkdir -p "$RUN" "$SCRATCH/tmp"
export TMPDIR="$SCRATCH/tmp/l5_oneway-${MODE}-${SLURM_JOB_ID:-local}"
mkdir -p "$TMPDIR"
cd "$RUN"

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_allow_ib=0
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export RUST_LOG="${RUST_LOG:-info}"
export PETSC_OPTIONS="-pc_factor_mat_solver_type mumps"
ulimit -n 65535 >/dev/null 2>&1 || true

echo "job ${SLURM_JOB_ID:-local} mode $MODE run dir $RUN  start $(date)"
echo "HEAD $(git -C "$REPO" rev-parse --short HEAD), dirty files: $(git -C "$REPO" status --porcelain -uno | wc -l)"

# Solid is the preCICE acceptor: start it first, then the BEM connector.
"$REPO/scripts/aeroenv.sh" aeroelast "$SOLID_YAML" --workdir "$RUN" > "$RUN/solid.log" 2>&1 &
SOLID_PID=$!
sleep 5
"$REPO/scripts/aeroenv.sh" aeroelast-bem-fsi "$FLUID_YAML" --workdir "$RUN" > "$RUN/fluid.log" 2>&1 &
FLUID_PID=$!

wait "$FLUID_PID"; FLUID_RC=$?
wait "$SOLID_PID"; SOLID_RC=$?
echo "fluid exit $FLUID_RC, solid exit $SOLID_RC, end $(date)"
tail -n 5 "$RUN/solid.log" "$RUN/fluid.log"
[ "$FLUID_RC" -eq 0 ] && [ "$SOLID_RC" -eq 0 ]
