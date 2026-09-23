#!/bin/bash
# submit_ch6_corotational.sh — one SLURM job per ch6 DLC prepared case (corotational only).
#
# Cases: dlc_1.1 (NTM rated), dlc_6.1 (EWM V50, parked), dlc_6.3 (EWM V1).
# Official WindIO blade (tests/IEA-15-240-RWT.yaml).
#
# Usage:
#   bash submit_ch6_corotational.sh
#   DLC_CASES="1.1 6.3" bash submit_ch6_corotational.sh
#   RESULTS_BASE=$SCRATCH/ch6_results_official bash submit_ch6_corotational.sh

set -euo pipefail

RESULTS_BASE="${RESULTS_BASE:-${SCRATCH}/ch6_results_official}"
read -ra DLC_CASES <<< "${DLC_CASES:-1.1 6.1 6.3}"

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASES_DIR="${SUBMIT_DIR}/ch6_prepared_cases"
RUN_SOLID="${SUBMIT_DIR}/_run_solid.sh"
RUN_FLUID="${SUBMIT_DIR}/_run_fluid.sh"

echo "Results base: ${RESULTS_BASE}"
echo "DLC cases:    ${DLC_CASES[*]}"

for f in "$RUN_SOLID" "$RUN_FLUID"; do
    [[ -x "$f" ]] || { echo "ERROR: not found or not executable: $f" >&2; exit 1; }
done

mkdir -p "${RESULTS_BASE}"

for dlc in "${DLC_CASES[@]}"; do
    solid_yaml="${CASES_DIR}/solid_corotational_dlc_${dlc}.yaml"
    fluid_yaml="${CASES_DIR}/fluid_dlc_${dlc}.yaml"
    for f in "$solid_yaml" "$fluid_yaml"; do
        [[ -f "$f" ]] || { echo "ERROR: missing YAML: $f" >&2; exit 1; }
    done

    case_dir="${RESULTS_BASE}/dlc_${dlc}"
    mkdir -p "${case_dir}/logs"

    job_name="ch6_${dlc//./_}_corot"
    sbatch \
        -p sequana_cpu \
        --ntasks=8 \
        -J "$job_name" \
        --time=4-00:00:00 \
        -e "${case_dir}/logs/slurm_%j.err" \
        -o "${case_dir}/logs/slurm_%j.out" \
        --mail-type=END,FAIL \
        --wrap="
set -e
cd \"${case_dir}\"
\"${RUN_SOLID}\" \"${solid_yaml}\" >> logs/solid.log 2>&1 &
SOLID_PID=\$!
sleep 2
\"${RUN_FLUID}\" \"${fluid_yaml}\" \"${case_dir}\" >> logs/fluid.log 2>&1 &
FLUID_PID=\$!
wait \$SOLID_PID
wait \$FLUID_PID
echo '[ch6 ${dlc}] done'
"
    echo "submitted ${job_name} -> ${case_dir}"
done
