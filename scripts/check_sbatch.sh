#!/usr/bin/env bash
# shellcheck shell=bash
#
# Submit scripts/check.sh to SLURM. Compute nodes have no git, so the revision and the file
# lists the gate needs (ruff's changed files, the mesh-size scan) are snapshotted here, on the
# login node, and handed to the job through CHECK_GIT_SNAPSHOT. The job tests the tree as it
# is when it starts; the snapshot records the tree at submission, so do not edit in between.
#
# Usage:
#   scripts/check_sbatch.sh full [--record] [-- extra sbatch options]
#   scripts/check_sbatch.sh quick [--record]
# Defaults: partition sequana_cpu, 240 min, 4 tasks (--ntasks only; no nodes, cpus-per-task or mem); log logs/check/<mode>-slurm-%j.log.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mode=full
record=""
sbatch_extra=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        quick | full) mode="$1" ;;
        --record) record=" --record" ;;
        --) shift; sbatch_extra=("$@"); break ;;
        -h | --help) sed -n '3,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "check_sbatch.sh: unknown argument '$1' (quick|full|--record|-- sbatch options)" >&2; exit 2 ;;
    esac
    shift
done

command -v git >/dev/null || { echo "check_sbatch.sh: git not on PATH; run this on the login node" >&2; exit 2; }

snap="$REPO_ROOT/logs/check/git-snapshot-$(date +%Y%m%dT%H%M%S)"
mkdir -p "$snap"
rev="$(git rev-parse --short HEAD)"
[ -z "$(git status --porcelain --untracked-files=no)" ] || rev="$rev+dirty"
echo "$rev" >"$snap/revision"
git ls-files --cached --others --exclude-standard tests tools src scripts >"$snap/tracked"
{
    git diff --name-only HEAD -- '*.py'
    git ls-files --others --exclude-standard -- '*.py'
} | sort -u >"$snap/changed_py"

sbatch --partition=sequana_cpu --time=240 --ntasks=4 \
    --job-name="check-$mode" --output="logs/check/$mode-slurm-%j.log" \
    --export="ALL,CHECK_GIT_SNAPSHOT=$snap,OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,RAYON_NUM_THREADS=4,PYTHONDONTWRITEBYTECODE=1" \
    ${sbatch_extra[@]+"${sbatch_extra[@]}"} \
    --wrap="scripts/check.sh $mode$record"
echo "check_sbatch.sh: snapshot ${snap#"$REPO_ROOT"/} (revision $rev)"
