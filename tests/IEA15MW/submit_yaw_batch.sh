#!/bin/bash
# submit_yaw_batch.sh — the whole yaw sweep as ONE Slurm job.
#
# One job per batch, with the cases running in parallel inside it: the cases do
# not parallelise individually (the BEM participant is serial, the solid uses 4
# OpenMP threads), so the parallelism has to come from running cases side by side.
# Each case gets its own CWD, which is what makes that safe — the coupling
# exchange directory is "." and precice-config.xml says so explicitly.
#
# Resources are requested with --ntasks only, never --mem.
#
# Usage:
#   bash submit_yaw_batch.sh                          # all 5 cases, corotational
#   SMOKE=1 bash submit_yaw_batch.sh                  # yaw_0 only, 2 s of physics
#   SMOKE=1 SMOKE_MAX_TIME=5 bash submit_yaw_batch.sh
#   SOLVER_TYPE=inertial bash submit_yaw_batch.sh
#   YAW_ANGLES="0" bash submit_yaw_batch.sh
#
# Output per case: ${RESULTS_BASE}/yaw_<angle>/{solid.log,fluid.log,logs/}

set -euo pipefail

SOLVER_TYPE="${SOLVER_TYPE:-corotational}"
SMOKE="${SMOKE:-0}"
SMOKE_MAX_TIME="${SMOKE_MAX_TIME:-2.0}"
read -ra YAW_ANGLES <<< "${YAW_ANGLES:-0 10 20 30 40}"
[[ "$SMOKE" == "1" ]] && YAW_ANGLES=("${SMOKE_YAW:-0}")

RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/frontiersin_results_${SOLVER_TYPE}_100s}"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SUBMIT_DIR}/frontiersin_2025_yaw"
RUN_SOLID="${SUBMIT_DIR}/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/_run_fluid.sh"
SMOKE_CASES_DIR="${SUBMIT_DIR}/yaw_smoke_cases"

# Cores per case: 4 OpenMP threads for the solid (the runner's default) + 1 for
# the serial BEM participant.  The whole batch is one job, so the request is
# cases x 5.
CORES_PER_CASE=5
NTASKS=$(( ${#YAW_ANGLES[@]} * CORES_PER_CASE ))

CONFIG="${SUBMIT_DIR}/precice-config.xml"
if [[ "$SMOKE" == "1" ]]; then
    # sits next to the real config so that "../precice-config-smoke.xml" and
    # the absolute path both resolve
    CONFIG="${SUBMIT_DIR}/precice-config-smoke.xml"
    sed -E "s/<max-time value=\"[0-9.]+\" \/>/<max-time value=\"${SMOKE_MAX_TIME}\" \/>/" \
        "${SUBMIT_DIR}/precice-config.xml" > "${CONFIG}"
    grep -q "max-time value=\"${SMOKE_MAX_TIME}\"" "${CONFIG}" || {
        echo "ERROR: could not set max-time in the smoke config" >&2; exit 1; }
fi

# ── Preflight ─────────────────────────────────────────────────────────────────
echo "Solver:       ${SOLVER_TYPE}"
echo "Results base: ${RESULTS_BASE}"
echo "Yaw angles:   ${YAW_ANGLES[*]} deg"
echo "Cases:        ${#YAW_ANGLES[@]}   ntasks: ${NTASKS} (${CORES_PER_CASE}/case)"
echo "preCICE cfg:  ${CONFIG}"
echo ""

for f in "$RUN_SOLID" "$RUN_FLUID"; do
    [[ -x "$f" ]] || { echo "ERROR: not found or not executable: $f" >&2; exit 1; }
done
for yaw in "${YAW_ANGLES[@]}"; do
    for f in "${CASES_DIR}/solid_${SOLVER_TYPE}_yaw_${yaw}.yaml" \
             "${CASES_DIR}/fluid_yaw_${yaw}.yaml"; do
        [[ -f "$f" ]] || { echo "ERROR: missing YAML: $f" >&2; exit 1; }
    done
done
echo "Preflight OK."

# ── Case directories, and a smoke copy of the YAMLs when asked ────────────────
declare -A SOLID_YAML FLUID_YAML CASE_DIR
for yaw in "${YAW_ANGLES[@]}"; do
    CASE_DIR[$yaw]="${RESULTS_BASE}/yaw_${yaw}"
    mkdir -p "${CASE_DIR[$yaw]}/logs"
    SOLID_YAML[$yaw]="${CASES_DIR}/solid_${SOLVER_TYPE}_yaw_${yaw}.yaml"
    FLUID_YAML[$yaw]="${CASES_DIR}/fluid_yaw_${yaw}.yaml"
    if [[ "$SMOKE" == "1" ]]; then
        # Point the case at the short preCICE config.  The copy has to live at the
        # SAME depth as frontiersin_2025_yaw/, because every other path in these
        # YAMLs is relative to the case file (yaml_file: ../../IEA-15-240-RWT.yaml,
        # airfoil_dir: ../../airfoils).  A flat copy resolves those one directory
        # too shallow and the run dies on a missing blade file — which is exactly
        # what the first smoke attempt did.
        mkdir -p "${SMOKE_CASES_DIR}"
        for pair in "solid:${SOLID_YAML[$yaw]}" "fluid:${FLUID_YAML[$yaw]}"; do
            kind="${pair%%:*}"; src="${pair#*:}"
            dst="${SMOKE_CASES_DIR}/$(basename "$src")"
            sed -E "s|config_file: .*precice-config[a-z-]*\.xml|config_file: ${CONFIG}|" \
                "$src" > "$dst"
            grep -q "config_file: ${CONFIG}" "$dst" || {
                echo "ERROR: could not repoint $src at the smoke config" >&2; exit 1; }
            if [[ "$kind" == "solid" ]]; then SOLID_YAML[$yaw]="$dst"; else FLUID_YAML[$yaw]="$dst"; fi
        done
    fi
done

# ── One job for the whole batch ───────────────────────────────────────────────
JOB_SCRIPT="$(mktemp "${RESULTS_BASE}/yaw_batch_XXXXXX.sh")"
{
    echo '#!/bin/bash'
    echo "#SBATCH --mail-type=END,FAIL"
    echo "#SBATCH --nodes=1"
    echo "#SBATCH --ntasks=${NTASKS}"
    echo "#SBATCH --time=24:00:00"
    echo "#SBATCH -e %j.err"
    echo "#SBATCH -J yaw_batch_${SOLVER_TYPE}"
    echo "#SBATCH -o %j.out"
    echo "#SBATCH -p sequana_cpu"
    echo
    echo 'echo "Job:    $SLURM_JOB_ID"'
    echo 'echo "Node:   $SLURM_NODELIST"'
    echo "echo \"Cases:  ${YAW_ANGLES[*]} (${#YAW_ANGLES[@]} in parallel, ${CORES_PER_CASE} cores each)\""
    echo 'echo "Start:  $(date)"'
    echo
    echo 'PIDS=(); NAMES=()'
    for yaw in "${YAW_ANGLES[@]}"; do
        echo "("
        echo "  cd '${CASE_DIR[$yaw]}'"
        echo "  '${RUN_SOLID}' '${SOLID_YAML[$yaw]}' > solid.log 2>&1"
        echo ") &"
        echo "PIDS+=(\$!); NAMES+=('solid:yaw_${yaw}')"
        echo "sleep 2"
        echo "("
        echo "  cd '${CASE_DIR[$yaw]}'"
        echo "  '${RUN_FLUID}' '${FLUID_YAML[$yaw]}' '${CASE_DIR[$yaw]}' > fluid.log 2>&1"
        echo ") &"
        echo "PIDS+=(\$!); NAMES+=('fluid:yaw_${yaw}')"
    done
    echo
    cat <<'BODY'
FAILED=0
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        echo "  OK   ${NAMES[$i]}"
    else
        echo "  FAIL ${NAMES[$i]} (exit=$?)" >&2
        FAILED=$((FAILED+1))
    fi
done
echo ""
echo "End: $(date)"
[[ $FAILED -eq 0 ]] && echo "all cases completed OK" || echo "$FAILED participant(s) failed"
exit $FAILED
BODY
} > "$JOB_SCRIPT"

JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
echo "Submitted batch job ${JOB_ID}  script: ${JOB_SCRIPT}"
echo "Monitor:  squeue -j ${JOB_ID}"
echo "Results:  ${RESULTS_BASE}/"
