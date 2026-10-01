#!/bin/bash
# DEPRECATED by tests/IEA15MW/submit_fsi_case_batch.sh (see case_lists/README.md):
# this script opens one Slurm job per case, which saturates the queue and, for the
# long cases, backfills far worse than one job per batch with the cases in parallel.
# It is kept for reference only; the case lists under case_lists/ are the
# reproducible way to launch these campaigns.
# submit_yaw_cases.sh — one SLURM job per yaw case (solid + fluid)
#
# Each job requests ~8 CPUs (solid uses 4 OMP threads + fluid uses 1).
# Jobs can run concurrently and move faster through the queue than a
# single 48-CPU job.
#
# Usage:
#   bash submit_yaw_cases.sh                        # all 5 cases, corotational
#   SOLVER_TYPE=inertial bash submit_yaw_cases.sh   # inertial solver
#   YAW_ANGLES="0 20" bash submit_yaw_cases.sh      # subset of angles
#   RESULTS_BASE=$SCRATCH/my_results bash submit_yaw_cases.sh
#
# Output layout (persistent):
#   $RESULTS_BASE/yaw_N/
#     corotational/ (or inertial/)  ← solid FEM output
#     fluid/                        ← BEM output
#     logs/
#       slurm_<JOB_ID>.out
#       slurm_<JOB_ID>.err
#       solid.log
#       fluid.log

set -euo pipefail

SOLVER_TYPE="${SOLVER_TYPE:-corotational}"
RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/frontiersin_results_${SOLVER_TYPE}_100s}"
read -ra YAW_ANGLES <<< "${YAW_ANGLES:-0 10 20 30 40}"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SUBMIT_DIR}/frontiersin_2025_yaw"
RUN_SOLID="${SUBMIT_DIR}/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/_run_fluid.sh"

# ── Preflight ──────────────────────────────────────────────────────────────────
echo "Solver:       ${SOLVER_TYPE}"
echo "Results base: ${RESULTS_BASE}"
echo "Yaw angles:   ${YAW_ANGLES[*]} deg"
echo ""

for f in "$RUN_SOLID" "$RUN_FLUID"; do
    [[ -x "$f" ]] || { echo "ERROR: not found or not executable: $f" >&2; exit 1; }
done

for yaw in "${YAW_ANGLES[@]}"; do
    for f in \
        "${CASES_DIR}/solid_${SOLVER_TYPE}_yaw_${yaw}.yaml" \
        "${CASES_DIR}/fluid_yaw_${yaw}.yaml"; do
        [[ -f "$f" ]] || { echo "ERROR: missing YAML: $f" >&2; exit 1; }
    done
done
echo "Preflight OK — all YAML files found."
echo ""

# ── Submit one job per case ────────────────────────────────────────────────────
SUBMITTED=()

for yaw in "${YAW_ANGLES[@]}"; do
    CASE_DIR="${RESULTS_BASE}/yaw_${yaw}"
    LOG_DIR="${CASE_DIR}/logs"
    SOLID_YAML="${CASES_DIR}/solid_${SOLVER_TYPE}_yaw_${yaw}.yaml"
    FLUID_YAML="${CASES_DIR}/fluid_yaw_${yaw}.yaml"

    mkdir -p "${LOG_DIR}"

    # sbatch reads the heredoc as the job script; --parsable returns only JOB_ID
    JOB_ID=$(sbatch --parsable << EOF
#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=24:00:00
#SBATCH -e %j.err
#SBATCH -J yaw_fsi_${yaw}deg_${SOLVER_TYPE}
#SBATCH -o %j.out
#SBATCH -p sequana_cpu

echo "Job:    \$SLURM_JOB_ID"
echo "Node:   \$SLURM_NODELIST"
echo "Case:   yaw=${yaw} deg  solver=${SOLVER_TYPE}"
echo "Output: ${CASE_DIR}"
echo "Start:  \$(date)"
echo ""

# solid is acceptor — start first and give it time to open the socket
( cd "${CASE_DIR}" && "${RUN_SOLID}" "${SOLID_YAML}" ) \
    >> "${LOG_DIR}/solid.log" 2>&1 &
SOLID_PID=\$!
echo "[yaw=${yaw}] solid PID=\${SOLID_PID}"

sleep 2

# fluid connects to solid
( cd "${CASE_DIR}" && "${RUN_FLUID}" "${FLUID_YAML}" "${CASE_DIR}" ) \
    >> "${LOG_DIR}/fluid.log" 2>&1 &
FLUID_PID=\$!
echo "[yaw=${yaw}] fluid PID=\${FLUID_PID}"

FAILED=0
wait "\$SOLID_PID" || { echo "[FAIL] solid (PID=\${SOLID_PID}) exit=\$?" >&2; FAILED=\$((FAILED+1)); }
wait "\$FLUID_PID" || { echo "[FAIL] fluid (PID=\${FLUID_PID}) exit=\$?" >&2; FAILED=\$((FAILED+1)); }

echo ""
echo "End: \$(date)"
[[ \$FAILED -eq 0 ]] && echo "✔ yaw=${yaw} completed OK" || echo "✖ \${FAILED} participant(s) failed"
exit \$FAILED
EOF
)

    echo "Submitted yaw=${yaw} deg → JOB_ID=${JOB_ID}  logs: ${LOG_DIR}/"
    SUBMITTED+=("${JOB_ID}:yaw_${yaw}")
done

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  ${#SUBMITTED[@]} job(s) submitted"
for entry in "${SUBMITTED[@]}"; do
    IFS=: read -r jid name <<< "$entry"
    echo "  ${jid}  ${name}"
done
echo ""
echo "Monitor:  squeue -u \$USER"
echo "Results:  ${RESULTS_BASE}/"
echo "════════════════════════════════════════════════════════════════"
