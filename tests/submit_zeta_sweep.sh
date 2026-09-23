#!/bin/bash
# submit_zeta_sweep.sh — one SLURM job per zeta case (solid + fluid)
#
# Zeta sweep: yaw=0 flexible FSI, 30 s, three Rayleigh damping ratios.
#   z_003   — zeta=0.03   (campaign baseline)
#   z_00048 — zeta=0.0048 (ElastoDyn BldFlDmp1, 0.48%)
#   z_000   — zeta=0.0    (no structural damping)
# Velocity feedback is inert in the corotational solver (solid never writes
# Velocity), so the fluid yaml uses velocity_data: null for all cases.
#
# Usage:
#   bash submit_zeta_sweep.sh
#   ZETA_CASES="z_003" bash submit_zeta_sweep.sh
#   RESULTS_BASE=$SCRATCH/zeta_sweep_results bash submit_zeta_sweep.sh

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SUBMIT_DIR}/zeta_sweep"
RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/zeta_sweep_results}"
RUN_SOLID="${SUBMIT_DIR}/IEA15MW/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/IEA15MW/_run_fluid.sh"

ZETA_CASES="${ZETA_CASES:-z_003 z_00048 z_000}"

for Z in ${ZETA_CASES}; do
    JOB_ID=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH -e %j.err
#SBATCH -J zeta_${Z}
#SBATCH -o %j.out
#SBATCH -p sequana_cpu

set -euo pipefail

Z="${Z}"
CASE_DIR="${RESULTS_BASE}/${Z}"
LOG_DIR="\${CASE_DIR}/logs"
mkdir -p "\${LOG_DIR}"

echo "Case: ${Z} (30 s, yaw=0)"
echo "Output: \${CASE_DIR}"

( cd "\${CASE_DIR}" && "${RUN_SOLID}" "${CASES_DIR}/${Z}/solid_corotational_yaw_0.yaml" ) \
    >> "\${LOG_DIR}/solid.log" 2>&1 &
SOLID_PID=\$!
sleep 2
( cd "\${CASE_DIR}" && "${RUN_FLUID}" "${CASES_DIR}/${Z}/fluid_yaw_0.yaml" "\${CASE_DIR}" ) \
    >> "\${LOG_DIR}/fluid.log" 2>&1 &
FLUID_PID=\$!

FAILED=0
wait "\${SOLID_PID}" || { echo "[FAIL] solid exit=\$?"; FAILED=\$((FAILED+1)); }
wait "\${FLUID_PID}" || { echo "[FAIL] fluid exit=\$?"; FAILED=\$((FAILED+1)); }

[[ \${FAILED} -eq 0 ]] && echo "ZETA ${Z} OK" || echo "ZETA ${Z} FAILED"
exit \${FAILED}
EOF
)
    echo "Submitted ${Z}: job ${JOB_ID}"
done
