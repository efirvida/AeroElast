#!/bin/bash
# submit_fsi_case_batch.sh — run a list of cases as ONE Slurm job.
#
# The same pattern as submit_yaw_batch.sh, generalised: the cases do not
# parallelise individually (the BEM participant is serial, the solid uses 4
# OpenMP threads), so one batch is one job and the parallelism comes from running
# cases side by side.  Resources are requested with --ntasks only, never --mem.
#
# Case list format, one case per line, `#` comments allowed:
#     <solid.yaml> <fluid.yaml> <case_dir>
#
# Usage:
#   bash submit_fsi_case_batch.sh cases.txt [-J jobname]
#   CORES_PER_CASE=5 MAX_CONCURRENT=9 bash submit_fsi_case_batch.sh cases.txt
#
# Cases run in waves of MAX_CONCURRENT: a sequana_cpu node has 48 cores, so 9
# cases at 5 cores each fit; a 15-case batch runs as 9 then 6 inside one job.
#
# -t/--time sets the wall limit.  Slurm's backfill uses the *requested* limit as
# the worst-case duration, so an oversized request is harder to place for no
# benefit -- but an undersized one silently truncates the run.  Size it from the
# measured per-window cost: at ~3.8 s per window a 100 s case needs ~11 h, while
# the parked stress-stiffened case runs at ~30 s per window and needs ~42 h.

set -euo pipefail

CASE_LIST="${1:?usage: submit_fsi_case_batch.sh <case-list> [-J jobname]}"
shift || true
JOB_NAME="fsi_batch"
TIME_LIMIT="24:00:00"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -J) JOB_NAME="$2"; shift 2 ;;
        -t) TIME_LIMIT="$2"; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

CORES_PER_CASE="${CORES_PER_CASE:-5}"
MAX_CONCURRENT="${MAX_CONCURRENT:-9}"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SOLID="${SUBMIT_DIR}/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/_run_fluid.sh"
[[ -x "$RUN_SOLID" && -x "$RUN_FLUID" ]] || { echo "ERROR: runners missing" >&2; exit 1; }
[[ -f "$CASE_LIST" ]] || { echo "ERROR: no case list at $CASE_LIST" >&2; exit 1; }

mapfile -t CASES < <(grep -vE '^\s*(#|$)' "$CASE_LIST")
N_CASES="${#CASES[@]}"
[[ "$N_CASES" -gt 0 ]] || { echo "ERROR: the case list is empty" >&2; exit 1; }
(( MAX_CONCURRENT > N_CASES )) && MAX_CONCURRENT="$N_CASES"
NTASKS=$(( MAX_CONCURRENT * CORES_PER_CASE ))

# ── Preflight: every file in every line must exist ────────────────────────────
for line in "${CASES[@]}"; do
    read -r solid fluid case_dir <<< "$line"
    for f in "$solid" "$fluid"; do
        [[ -f "$f" ]] || { echo "ERROR: missing case file: $f" >&2; exit 1; }
    done
    mkdir -p "$case_dir"
done
echo "Cases:        ${N_CASES}"
echo "Concurrency:  ${MAX_CONCURRENT} (waves of ${MAX_CONCURRENT}), ${CORES_PER_CASE} cores/case"
echo "ntasks:       ${NTASKS}"
echo "Job name:     ${JOB_NAME}   time limit: ${TIME_LIMIT}"

JOB_SCRIPT="$(mktemp "/tmp/${JOB_NAME}_XXXXXX.sh")"
{
    echo '#!/bin/bash'
    echo "#SBATCH --mail-type=END,FAIL"
    echo "#SBATCH --nodes=1"
    echo "#SBATCH --ntasks=${NTASKS}"
    echo "#SBATCH --time=${TIME_LIMIT}"
    echo "#SBATCH -e %j.err"
    echo "#SBATCH -J ${JOB_NAME}"
    echo "#SBATCH -o %j.out"
    echo "#SBATCH -p sequana_cpu"
    echo
    echo 'echo "Job:   $SLURM_JOB_ID  node: $SLURM_NODELIST"'
    echo 'echo "Start: $(date)"'
    echo "MAX_CONCURRENT=${MAX_CONCURRENT}"
    echo "SOLID='${RUN_SOLID}'; FLUID='${RUN_FLUID}'"
    echo 'FAILED=0; DONE=0; RUNNING=0'
    echo 'pids=(); names=(); dirs=()'
    echo 'reap() {'
    echo '  local i'
    echo '  for i in "${!pids[@]}"; do'
    echo '    if [[ -n "${pids[$i]}" ]] && ! kill -0 "${pids[$i]}" 2>/dev/null; then'
    echo '      if wait "${pids[$i]}"; then echo "  OK   ${names[$i]}"; else echo "  FAIL ${names[$i]}" >&2; FAILED=$((FAILED+1)); fi'
    echo '      pids[$i]=""; RUNNING=$((RUNNING-1))'
    echo '    fi'
    echo '  done'
    echo '}'
    echo 'launch() {'
    echo '  # The two participants are paired: when either exits, the other is killed.'
    echo '  # A dead solid leaves the BEM waiting on preCICE forever, and the job'
    echo '  # would then hold the node until the wall limit -- one failed case must'
    echo '  # not cost 24 h of a queue slot.' 
    echo '  local solid="$1" fluid="$2" dir="$3" tag="$4"'
    echo '  mkdir -p "$dir"'
    echo '  ('
    echo '    cd "$dir" || exit 1'
    echo '    "$SOLID" "$solid" > solid.log 2>&1 &'
    echo '    local sp=$!'
    echo '    ( sleep 2; "$FLUID" "$fluid" "$dir" > fluid.log 2>&1 ) &'
    echo '    local fp=$!'
    echo '    wait -n "$sp" "$fp"; local rc=$?'
    echo '    kill "$sp" "$fp" 2>/dev/null'
    echo '    wait 2>/dev/null'
    echo '    exit $rc'
    echo '  ) &'
    echo '  pids+=($!); names+=("case:$tag"); dirs+=("$dir"); RUNNING=$((RUNNING+1))'
    echo '}'
    echo
    echo 'while IFS= read -r line; do'
    echo '  [[ -z "${line// }" || "$line" == \#* ]] && continue'
    echo '  read -r solid fluid dir <<< "$line"'
    echo '  while (( RUNNING >= MAX_CONCURRENT * 2 )); do reap; sleep 5; done'
    echo '  launch "$solid" "$fluid" "$dir" "$(basename "$dir")"'
    echo 'done <<'"'"'CASES'"'"''
    printf '%s\n' "${CASES[@]}"
    echo 'CASES'
    echo 'while (( RUNNING > 0 )); do reap; sleep 5; done'
    echo 'echo "End: $(date)"'
    echo '[[ $FAILED -eq 0 ]] && echo "all cases completed OK" || echo "$FAILED participant(s) failed"'
    echo 'exit $FAILED'
} > "$JOB_SCRIPT"

JOB_ID=$(sbatch --parsable "$JOB_SCRIPT")
echo "Submitted job ${JOB_ID}  script: ${JOB_SCRIPT}"
echo "Monitor: squeue -j ${JOB_ID}"
