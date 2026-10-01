#!/bin/bash
# DEPRECATED by tests/IEA15MW/submit_fsi_case_batch.sh (see case_lists/README.md):
# this script opens one Slurm job per case, which saturates the queue and, for the
# long cases, backfills far worse than one job per batch with the cases in parallel.
# It is kept for reference only; the case lists under case_lists/ are the
# reproducible way to launch these campaigns.
# submit_parked_v50.sh — parked DLC 6.x-analog case (V50, yaw=8°, pitch=90°).
# Official WindIO blade. Inputs preserved from the old bem_90_50_S case.
#
# Usage:
#   bash submit_parked_v50.sh
#   RESULTS_BASE=$SCRATCH/parked_v50_official bash submit_parked_v50.sh

set -euo pipefail

RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/parked_v50_official}"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="${SUBMIT_DIR}/parked_v50"
RUN_SOLID="${SUBMIT_DIR}/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/_run_fluid.sh"

for f in "$RUN_SOLID" "$RUN_FLUID" "$CASE_DIR/solid_v50.yaml" "$CASE_DIR/fluid_v50.yaml"; do
    [[ -e "$f" ]] || { echo "ERROR: missing: $f" >&2; exit 1; }
done

mkdir -p "${RESULTS_BASE}/logs"

sbatch \
    -p sequana_cpu \
    --ntasks=8 \
    -J parked_v50_official \
    --time=4-00:00:00 \
    -e "${RESULTS_BASE}/logs/slurm_%j.err" \
    -o "${RESULTS_BASE}/logs/slurm_%j.out" \
    --mail-type=END,FAIL \
    --wrap="
set -e
cd \"${RESULTS_BASE}\"
\"${RUN_SOLID}\" \"${CASE_DIR}/solid_v50.yaml\" >> logs/solid.log 2>&1 &
SOLID_PID=\$!
sleep 2
\"${RUN_FLUID}\" \"${CASE_DIR}/fluid_v50.yaml\" \"${RESULTS_BASE}\" >> logs/fluid.log 2>&1 &
FLUID_PID=\$!
wait \$SOLID_PID
wait \$FLUID_PID
echo '[parked_v50] done'
"
echo "submitted parked_v50 -> ${RESULTS_BASE}"
