"""Inspect reduced-order VLM sectional diagnostics for one operating point."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from aeroelast.solvers.aero.validation import (
    RotorOperatingPoint,
    compute_vlm_sectional_diagnostics,
)


def _resolve_default_blade_file(repo_root: Path) -> Path:
    candidates = [
        repo_root / "tests" / "IEA-15-240-RWT.yaml",
        repo_root / "examples" / "reference_turbines" / "yamls" / "IEA-15-240-RWT.yaml",
        repo_root.parent / "simulations" / "blade" / "solid" / "IEA-15-240-RWT.yaml",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "IEA-15-240-RWT.yaml not found. Use --blade-file to point at the WindIO blade definition."
    )


def _rows_from_diagnostics(diagnostics: object) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    count = len(diagnostics.radii)
    for i in range(count):
        rows.append({
            "radius_m": float(diagnostics.radii[i]),
            "alpha_incident_deg": float(diagnostics.alpha_incident_deg[i]),
            "alpha_total_deg": float(diagnostics.alpha_total_deg[i]),
            "cl_polar": float(diagnostics.cl_polar[i]),
            "cd_polar": float(diagnostics.cd_polar[i]),
            "cl_vlm_trailing": float(diagnostics.cl_vlm_trailing[i]),
            "gamma_leading": float(diagnostics.gamma_leading[i]),
            "gamma_trailing": float(diagnostics.gamma_trailing[i]),
            "strip_base_thrust_N": float(diagnostics.strip_base_thrust[i]),
            "strip_base_torque_Nm": float(diagnostics.strip_base_torque[i]),
        })
    return rows


def _print_rows(rows: list[dict[str, float]]) -> None:
    print(
        "radius[m]  alpha_inc[deg]  alpha_tot[deg]  cl_polar  cl_vlm_te  "
        "gamma_le  gamma_te  dT_base[N]  dQ_base[Nm]"
    )
    for row in rows:
        print(
            f"{row['radius_m']:8.3f}"
            f"{row['alpha_incident_deg']:16.3f}"
            f"{row['alpha_total_deg']:15.3f}"
            f"{row['cl_polar']:10.3f}"
            f"{row['cl_vlm_trailing']:11.3f}"
            f"{row['gamma_leading']:10.3f}"
            f"{row['gamma_trailing']:10.3f}"
            f"{row['strip_base_thrust_N']:12.3f}"
            f"{row['strip_base_torque_Nm']:13.3f}"
        )


def _write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blade-file", type=Path, help="Path to the WindIO blade YAML file.")
    parser.add_argument("--name", default="sectional-diagnostic")
    parser.add_argument("--wind-speed", type=float, default=10.659)
    parser.add_argument("--omega-rpm", type=float, default=7.559)
    parser.add_argument("--pitch-deg", type=float, default=0.0)
    parser.add_argument("--n-spanwise-panels", type=int, default=12)
    parser.add_argument("--n-chordwise-panels", type=int, default=2)
    parser.add_argument("--csv", type=Path, help="Optional CSV output path.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    blade_file = (
        args.blade_file.resolve() if args.blade_file else _resolve_default_blade_file(repo_root)
    )

    operating_point = RotorOperatingPoint(
        name=args.name,
        wind_speed=args.wind_speed,
        omega_rpm=args.omega_rpm,
        pitch_deg=args.pitch_deg,
    )

    diagnostics = compute_vlm_sectional_diagnostics(
        blade_file,
        operating_point,
        n_spanwise_panels=args.n_spanwise_panels,
        n_chordwise_panels=args.n_chordwise_panels,
    )

    rows = _rows_from_diagnostics(diagnostics)
    _print_rows(rows)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        _write_csv(args.csv, rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
