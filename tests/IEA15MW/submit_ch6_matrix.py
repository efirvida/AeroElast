#!/usr/bin/env python3
"""V-06 full DLC 1.1 matrix launcher: the 12 wind bins x seeds + extra points.

Generates the per-case YAMLs (the wind / pitch / omega injected from the IEA
15 MW operating schedule) and submits one SLURM job per case with the 96-h
wall (the maximum for the queue).

Operating schedule (the simplified IEA 15 MW controller):
  * below rated (wind <= 10.59 m/s): TSR-9 omega, capped at the rated
    omega = 0.7872 rad/s (7.52 rpm); pitch = 0
  * above rated: the rated omega, the pitch ramping ~1.8 deg per m/s

NOTE: the fluid participant is the BEM with the uniform inflow, so the
"seeds" produce identical runs (the NTM turbulence is not modelled — the
turbulent wind fields are a separate future work item).  The launcher
supports them for the matrix completeness; drop them with SEEDS="1" to
save the compute.

Usage:
  SEEDS="1" python tests/IEA15MW/submit_ch6_matrix.py     # 12 runs
  python tests/IEA15MW/submit_ch6_matrix.py               # 72 + 3 extra
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CASES_DIR = REPO / "tests" / "IEA15MW" / "ch6_prepared_cases"
GEN_DIR = CASES_DIR  # the generated yamls live at the same level so the
                     # relative asset paths (../../IEA-15-240-RWT.yaml etc.)
                     # keep resolving exactly like the templates
RESULTS_BASE = Path(os.environ.get("RESULTS_BASE",
                                   "/scratch/leahk/eduardo.donestevez/v006_matrix_results"))

WIND_SPEEDS = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]
EXTRA_POINTS = [8, 10, 12]  # the finer FSI operating points around the rated
SEEDS = [int(s) for s in os.environ.get("SEEDS", "1 2 3 4 5 6").split()]

RATED_WIND = 10.59
RATED_OMEGA = 0.7872360607268416  # rad/s (7.52 rpm)
R = 120.97  # m
TSR_OPT = 9.0
PITCH_SLOPE = 1.8  # deg per m/s above the rated


def schedule(wind: float) -> tuple[float, float]:
    omega = min(TSR_OPT * wind / R, RATED_OMEGA)
    pitch = max(0.0, (wind - RATED_WIND) * PITCH_SLOPE)
    return omega, pitch


def make_case(case_id: str, wind: float, seed: int) -> tuple[Path, Path]:
    omega, pitch = schedule(wind)
    solid_src = CASES_DIR / "solid_corotational_dlc_1.1.yaml"
    fluid_src = CASES_DIR / "fluid_dlc_1.1.yaml"
    solid_dst = GEN_DIR / f"gen_{case_id}_solid.yaml"
    fluid_dst = GEN_DIR / f"gen_{case_id}_fluid.yaml"

    s = solid_src.read_text()
    s = s.replace("    omega: 0.7872360607268416", f"    omega: {omega}")
    s = s.replace("  flow_velocity: 10.8", f"  flow_velocity: {wind}")
    solid_dst.write_text(s)

    f = fluid_src.read_text()
    f = f.replace("  wind_speed: 10.8", f"  wind_speed: {wind}")
    f = f.replace("  pitch: 0.0", f"  pitch: {pitch}")
    f = f.replace("    initial: 0.7872360607268416", f"    initial: {omega}")
    fluid_dst.write_text(f)
    return solid_dst, fluid_dst


def submit(case_id: str, solid_yaml: Path, fluid_yaml: Path) -> str:
    case_dir = RESULTS_BASE / case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "logs").mkdir(exist_ok=True)
    marker = case_dir / ".submitted"
    if marker.exists():
        return "SKIPPED"
    run_solid = REPO / "tests" / "IEA15MW" / "_run_solid.sh"
    run_fluid = REPO / "tests" / "IEA15MW" / "_run_fluid.sh"
    wrap = (
        f'set -e\ncd "{case_dir}"\n'
        f'"{run_solid}" "{solid_yaml}" >> logs/solid.log 2>&1 &\n'
        f"SOLID_PID=$!\nsleep 2\n"
        f'"{run_fluid}" "{fluid_yaml}" "{case_dir}" >> logs/fluid.log 2>&1 &\n'
        f"FLUID_PID=$!\nwait $SOLID_PID\nwait $FLUID_PID\n"
        f"echo '[v006 {case_id}] done'"
    )
    # pace against the queue's MaxSubmitPU: retry on the submit limit
    for attempt in range(120):
        res = subprocess.run(
            ["sbatch", "--parsable", "-p", "sequana_cpu", "--ntasks=8",
             "-J", f"v006_{case_id}", "--time=4-00:00:00",
             "-e", f"{case_dir}/logs/slurm_%j.err",
             "-o", f"{case_dir}/logs/slurm_%j.out",
             "--mail-type=END,FAIL", "--wrap", wrap],
            capture_output=True, text=True)
        if res.returncode == 0:
            marker.write_text(res.stdout.strip())
            return res.stdout.strip()
        if "MaxSubmit" not in res.stderr and "MaxJobs" not in res.stderr:
            raise RuntimeError(f"sbatch failed for {case_id}: {res.stderr}")
        time.sleep(60)
    raise RuntimeError(f"submit limit not released after 2 h for {case_id}")


def main() -> int:
    GEN_DIR.mkdir(exist_ok=True)
    submitted = []

    for wind in WIND_SPEEDS:
        for seed in SEEDS:
            case_id = f"v{wind:g}_s{seed}"
            solid, fluid = make_case(case_id, wind, seed)
            jid = submit(case_id, solid, fluid)
            submitted.append((case_id, jid))
            print(f"{case_id} -> {jid}", flush=True)

    for wind in EXTRA_POINTS:
        case_id = f"op_v{wind:g}"
        solid, fluid = make_case(case_id, wind, 1)
        jid = submit(case_id, solid, fluid)
        submitted.append((case_id, jid))
        print(f"{case_id} -> {jid}", flush=True)

    print(f"\ntotal: {len(submitted)} jobs -> {RESULTS_BASE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
