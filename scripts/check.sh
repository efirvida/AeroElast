#!/usr/bin/env bash
# shellcheck shell=bash
#
# Repo-local verification gate for fem-shell.
#
# The full pytest suite takes 20-26 minutes and ends with 190k-270k warnings, so it is a
# measurement run, not a smoke test. This gate exists so that "did I break something?" is
# a ~15 second question, and so that the documented-red tests stop being re-diagnosed from
# scratch every session (see scripts/known-red-<mode>.txt; the repo carries four
# documented-red tests, and the quick selection currently includes one of them).
#
# Usage:
#   scripts/check.sh quick [--record]   # Ruff on changed files + fast physical subset
#   scripts/check.sh full  [--record]   # the documented suite, quiet, log kept
#   scripts/check.sh --list             # show the selections and the known-red registry
#
# Exit codes: 0 no new failures, 1 a failure outside the mode's registry, 2 the environment
# or pytest could not start. A known-red failure alone is exit 0 — that is the point of the
# registry — but it is always printed.
#
# This is a smoke gate. It does not replace the physical validation anchors (S-2, S-4,
# S-5, S-7) or a full-suite measurement before a merge.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# One registry per mode: a full-suite measurement must not clobber the quick registry.
KNOWN_RED=""
LOG_DIR="$REPO_ROOT/logs/check"
mkdir -p "$LOG_DIR"

# Fast physical subset: element/DOF conventions, laminate coupling, omega rebuild,
# section contour, rotor aero geometry, the box torsion benchmark, and the corotational
# tangent check. Measured at ~13 s together on the 2026-10-09 tree.
QUICK_TESTS=(
    tests/test_mixed_mesh_convention.py
    tests/test_omega_rebuild_predicate.py
    tests/test_laminate_bend_twist.py
    tests/test_laminate_shear_coupling.py
    tests/test_section_contour.py
    tests/test_rotor_aero_geometry.py
    tests/test_box_torsion_bending_benchmark.py
    tests/test_corotational_large_rotation_validation.py
)

# Machine-readable pytest: the project addopts (`-v --tb=short`) triples the output of a
# passing run and prints the per-test list nobody acts on. `--disable-warnings` drops the
# summary block; the warning count still appears in the final line.
PYTEST_OPTS=(-o addopts="" -q --tb=line --disable-warnings -p no:cacheprovider)

mode=""
record=0
for arg in "$@"; do
    case "$arg" in
        quick | full) mode="$arg" ;;
        --record) record=1 ;;
        --list) mode="list" ;;
        -h | --help)
            sed -n '3,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "check.sh: unknown argument '$arg' (quick|full|--record|--list)" >&2
            exit 2
            ;;
    esac
done
mode="${mode:-quick}"
KNOWN_RED="$REPO_ROOT/scripts/known-red-$mode.txt"

read_known_red() {
    [ -f "$KNOWN_RED" ] || return 0
    sed -e 's/#.*//' -e 's/[[:space:]]*$//' "$KNOWN_RED" | grep -v '^$' || true
}

if [ "$mode" = "list" ]; then
    echo "quick selection (${#QUICK_TESTS[@]} files):"
    printf '  %s\n' "${QUICK_TESTS[@]}"
    echo
    for candidate in "$REPO_ROOT"/scripts/known-red-*.txt; do
        [ -f "$candidate" ] || continue
        echo "$(basename "$candidate"):"
        KNOWN_RED="$candidate" read_known_red | sed 's/^/  /'
    done
    exit 0
fi

log="$LOG_DIR/$mode-$(date +%Y%m%dT%H%M%S).log"

# ---- environment ---------------------------------------------------------------------
# shellcheck source=scripts/aeroenv.sh
if ! source "$REPO_ROOT/scripts/aeroenv.sh"; then
    echo "check: FAIL — aeroelast runtime could not be assembled; run scripts/aeroenv.sh --check" >&2
    exit 2
fi

# ---- Ruff: only the files this change touched -----------------------------------------
# A repo-wide `ruff check` reports 313 pre-existing errors (161 kB of output). pi-lens
# already lints edited files, so the gate checks the diff and says nothing when it is clean.
changed_py=()
while IFS= read -r path; do
    [ -n "$path" ] && [ -f "$path" ] && changed_py+=("$path")
done < <(
    {
        git diff --name-only HEAD -- '*.py' 2>/dev/null
        git ls-files --others --exclude-standard -- '*.py' 2>/dev/null
    } | sort -u
)

ruff_status="clean"
if [ "${#changed_py[@]}" -gt 0 ]; then
    if ! ruff check --output-format=concise "${changed_py[@]}" >"$log.ruff" 2>&1; then
        ruff_status="issues in ${#changed_py[@]} changed file(s)"
    fi
else
    ruff check --statistics >"$log.ruff" 2>&1
    ruff_status="no changed .py files (repo total: $(grep -c $'\t' "$log.ruff" 2>/dev/null || echo 0) rule buckets)"
fi

# ---- pytest ---------------------------------------------------------------------------
if [ "$mode" = "quick" ]; then
    selection=("${QUICK_TESTS[@]}")
else
    selection=(tests/)
fi

python -m pytest "${PYTEST_OPTS[@]}" "${selection[@]}" >"$log" 2>&1
pytest_rc=$?

summary="$(grep -E '^[0-9]+ (passed|failed|error)|=+ .* (passed|failed|error).* =+$' "$log" | tail -1)"
[ -n "$summary" ] || summary="$(tail -1 "$log")"

failed_ids=()
while IFS= read -r line; do
    [ -n "$line" ] && failed_ids+=("${line#FAILED }")
done < <(grep -E '^FAILED ' "$log" || true)

# ---- record ---------------------------------------------------------------------------
# Record before classifying: a recorded run must not also report its own entries as new.
if [ "$record" = "1" ]; then
    {
        echo "# Known-red registry for scripts/check.sh (mode: $mode)."
        echo "# Re-derived automatically by: scripts/check.sh $mode --record"
        echo "# Recorded: $(date -Iseconds) at revision $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
        echo "# Invalidated by: any change under crates/, a maturin rebuild, or a new failure"
        echo "# that appears here outside this list. Never hand-edit to make a run green."
        printf '%s\n' ${failed_ids[@]+"${failed_ids[@]}"}
    } >"$KNOWN_RED"
    echo "check[$mode]: recorded ${#failed_ids[@]} known-red entry(ies) in ${KNOWN_RED#"$REPO_ROOT"/}"
fi

# ---- classify -------------------------------------------------------------------------
known="$(read_known_red)"
new_failures=()
known_failures=()
for id in ${failed_ids[@]+"${failed_ids[@]}"}; do
    if printf '%s\n' "$known" | grep -qxF -- "$id"; then
        known_failures+=("$id")
    else
        new_failures+=("$id")
    fi
done

# A registry entry that now passes is stale knowledge, and stale knowledge is the thing
# this gate is supposed to remove.
healed=()
while IFS= read -r id; do
    [ -n "$id" ] || continue
    printf '%s\n' ${failed_ids[@]+"${failed_ids[@]}"} | grep -qxF -- "$id" || healed+=("$id")
done <<<"$known"

# ---- report ---------------------------------------------------------------------------
echo "check[$mode]: ruff: $ruff_status"
echo "check[$mode]: pytest: ${summary:-no summary}"
echo "check[$mode]: log: ${log#"$REPO_ROOT"/}"
for id in ${known_failures[@]+"${known_failures[@]}"}; do
    echo "  KNOWN-RED  $id"
done
for id in ${new_failures[@]+"${new_failures[@]}"}; do
    echo "  NEW-FAIL   $id"
done
for id in ${healed[@]+"${healed[@]}"}; do
    echo "  now-passes $id   (re-record: scripts/check.sh $mode --record)"
done

if [ "${#new_failures[@]}" -gt 0 ]; then
    echo "check[$mode]: FAIL — ${#new_failures[@]} failure(s) outside ${KNOWN_RED#"$REPO_ROOT"/}"
    exit 1
fi
if [ "$pytest_rc" -gt 1 ] || [ -z "$summary" ]; then
    echo "check[$mode]: FAIL — pytest did not complete (rc=$pytest_rc); see ${log#"$REPO_ROOT"/}"
    exit 2
fi
echo "check[$mode]: OK"
exit 0
