#!/usr/bin/env bash
# shellcheck shell=bash
#
# Runtime bootstrap for the aeroelast Python + Rust hybrid.
#
# The Rust extension is built with GCC 14, so importing `aeroelast` (or `_aeroelast`)
# against the system libstdc++ fails with:
#
#     ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
#
# Shell state does not survive between agent tool calls, so every command that touches
# the package needs this bootstrap — not just the first one.
#
# Usage:
#   source scripts/aeroenv.sh                 # export into the current shell
#   scripts/aeroenv.sh <command> [args...]    # run <command> with a valid environment
#   scripts/aeroenv.sh --check                # verify the environment and print one line
#
# Environment overrides:
#   VENV_DIR            virtualenv root (default: <repo>/.venv)
#   AEROENV_GCC_MODULE  environment module to load (default: gcc/14.2.0_sequana)
#   AEROENV_VERBOSE=1   print where the runtime was found (stderr)
#
# Exit codes: 0 ok, 2 the environment could not be assembled, otherwise the wrapped
# command's own exit code.

_aeroenv_ok=1

_aeroenv_find_repo_root() {
    local here="${BASH_SOURCE[0]}"
    while [ -L "$here" ]; do
        local target
        target="$(readlink "$here")"
        case "$target" in
            /*) here="$target" ;;
            *) here="$(dirname "$here")/$target" ;;
        esac
    done
    (cd "$(dirname "$here")/.." && pwd)
}

_aeroenv_has_symbol() {
    [ -e "$1" ] || return 1
    # No pipe on purpose: under the caller's `set -o pipefail`, `strings | grep -q` fails
    # on grep's early exit (SIGPIPE) even when the symbol is present.
    local symbols
    symbols="$(strings "$1" 2>/dev/null || true)"
    case "$symbols" in
        *CXXABI_1.3.15*) return 0 ;;
    esac
    return 1
}

# Drop empty and repeated entries: the container's LD_LIBRARY_PATH already carries the
# venv paths, and `module load` relocates the GCC tree, so naive concatenation grows.
_aeroenv_dedupe_path() {
    local result="" rest="$1" entry
    local -A seen=()
    while [ -n "$rest" ]; do
        case "$rest" in
            *:*) entry="${rest%%:*}"; rest="${rest#*:}" ;;
            *) entry="$rest"; rest="" ;;
        esac
        [ -n "$entry" ] || continue
        [ -n "${seen[$entry]:-}" ] && continue
        seen[$entry]=1
        result="${result:+$result:}$entry"
    done
    printf '%s' "$result"
}

_aeroenv_setup() {
    local repo_root="$1"
    local venv="${VENV_DIR:-$repo_root/.venv}"
    local module_name="${AEROENV_GCC_MODULE:-gcc/14.2.0_sequana}"

    if [ ! -x "$venv/bin/python" ]; then
        echo "aeroenv: no python at $venv/bin/python (set VENV_DIR)" >&2
        return 2
    fi

    # 1. GCC 14 runtime. Prefer the environment module; fall back to the installed tree.
    if type -t module >/dev/null 2>&1; then
        module load "$module_name" >/dev/null 2>&1 || true
    fi

    local gcc_lib=""
    local candidate
    candidate="$(gcc -print-file-name=libstdc++.so.6 2>/dev/null)"
    if [ -n "$candidate" ] && _aeroenv_has_symbol "$candidate"; then
        gcc_lib="$(dirname "$candidate")"
    else
        for candidate in \
            "/scratch/app_sequana/gcc/14.2.0/lib64" \
            "/scratch/app_sequana/gcc/14.2.0/lib" \
            "/petrobr/app_sequana/gcc/14.2.0/lib64" \
            "/petrobr/app_sequana/gcc/14.2.0/lib"; do
            if _aeroenv_has_symbol "$candidate/libstdc++.so.6"; then
                gcc_lib="$candidate"
                break
            fi
        done
    fi

    if [ -z "$gcc_lib" ]; then
        echo "aeroenv: no libstdc++ exporting CXXABI_1.3.15; tried \`module load $module_name\`" >&2
        echo "aeroenv: and the known GCC 14 prefixes. Set AEROENV_GCC_MODULE or fix the module path." >&2
        return 2
    fi
    gcc_lib="$(cd "$gcc_lib" 2>/dev/null && pwd -P || printf '%s' "$gcc_lib")"

    # 2. LD_LIBRARY_PATH, rebuilt from a saved base so repeated sourcing is idempotent.
    local lib_path="$gcc_lib:$venv/lib:$venv/lib64"
    [ -d /usr/lib64/psm2-compat ] && lib_path="$lib_path:/usr/lib64/psm2-compat"
    if [ -n "${AEROENV_BASE_LD_LIBRARY_PATH:-}" ]; then
        lib_path="$lib_path:$AEROENV_BASE_LD_LIBRARY_PATH"
    elif [ -n "${LD_LIBRARY_PATH:-}" ]; then
        AEROENV_BASE_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
        lib_path="$lib_path:$LD_LIBRARY_PATH"
    fi

    lib_path="$(_aeroenv_dedupe_path "$lib_path")"

    # 3. venv first on PATH, and VIRTUAL_ENV for maturin.
    case ":$PATH:" in
        *":$venv/bin:"*) ;;
        *) PATH="$venv/bin:$PATH" ;;
    esac
    export PATH
    export VIRTUAL_ENV="$venv"
    export LD_LIBRARY_PATH="$lib_path"
    export AEROENV_GCC_LIB="$gcc_lib"

    if [ "${AEROENV_VERBOSE:-0}" = "1" ]; then
        echo "aeroenv: gcc14=$gcc_lib venv=$venv" >&2
    fi
    return 0
}

_aeroenv_repo_root="$(_aeroenv_find_repo_root)"

if _aeroenv_setup "$_aeroenv_repo_root"; then
    _aeroenv_ok=0
fi

# Sourced: leave the environment in place and return.
if [ "${BASH_SOURCE[0]}" != "${0}" ]; then
    # shellcheck disable=SC2317  # `return` legitimately fails when executed, not sourced
    return "$_aeroenv_ok" 2>/dev/null || exit "$_aeroenv_ok"
fi

if [ "$_aeroenv_ok" != "0" ]; then
    exit 2
fi

if [ "$#" -eq 0 ] || [ "$1" = "--check" ]; then
    echo "aeroenv: ok ($(python -c 'import sys; print(sys.executable)') · gcc14=$AEROENV_GCC_LIB)"
    exit 0
fi

# An absolute path is required: `exec` on the bare name would re-search PATH, which is
# fine, but using the venv interpreter directly is what keeps `python` honest here.
case "$1" in
    python | python3) exec "$VIRTUAL_ENV/bin/python" "${@:2}" ;;
    *) exec "$@" ;;
esac
