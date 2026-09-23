#!/bin/bash
# submit_exp_feedback_off.sh — one job per feedback-isolation case (solid + fluid)
#
# Isolates the two geometric feedback terms of the deformed BEM rebuild:
#   twist_off   — deformed_twist: false (twist_def = ref_twist)
#   radius_off  — deformed_radius: false (r_def = ref_r)
# Both: yaw=0 flexible FSI, 30 s, zeta=0.03, velocity feedback already inert.
#
# Usage:
#   bash submit_exp_feedback_off.sh
#   EXP_CASES="twist_off" bash submit_exp_feedback_off.sh

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SUBMIT_DIR}/exp_feedback_off"
RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/exp_feedback_off_results}"
RUN_SOLID="${SUBMIT_DIR}/IEA15MW/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/IEA15MW/_run_fluid.sh"

EXP_CASES="${EXP_CASES:-twist_off radius_off}"

for E in ${EXP_CASES}; do
    JOB_ID=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH -e %j.err
#SBATCH -J fb_${E}
#SBATCH -o %j.out
#SBATCH -p sequana_cpu

set -euo pipefail

E="${E}"
CASE_DIR="${RESULTS_BASE}/${E}"
LOG_DIR="\${CASE_DIR}/logs"
mkdir -p "\${LOG_DIR}"

echo "Case: \${E} (30 s, yaw=0)"
echo "Output: \${CASE_DIR}"

( cd "\${CASE_DIR}" && "${RUN_SOLID}" "${CASES_DIR}/${E}/solid_corotational_yaw_0.yaml" ) \
    >> "\${LOG_DIR}/solid.log" 2>&1 &
SOLID_PID=\$!
sleep 2
( cd "\${CASE_DIR}" && "${RUN_FLUID}" "${CASES_DIR}/${E}/fluid_yaw_0.yaml" "\${CASE_DIR}" ) \
    >> "\${LOG_DIR}/fluid.log" 2>&1 &
FLUID_PID=\$!

FAILED=0
wait "\${SOLID_PID}" || { echo "[FAIL] solid exit=\$?"; FAILED=\$((FAILED+1)); }
wait "\${FLUID_PID}" || { echo "[FAIL] fluid exit=\$?"; FAILED=\$((FAILED+1)); }

[[ \${FAILED} -eq 0 ]] && echo "FB \${E} OK" || echo "FB \${E} FAILED"
exit \${FAILED}
EOF
)
    echo "Submitted ${E}: job ${JOB_ID}"
done
