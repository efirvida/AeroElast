"""A7 — Cross-baseline BEM-only torque/power at V-01 and V-04 operating points.

Resolves the rigid-baseline ambiguity in the corotational report: the rigid baseline
reported in §3.7.1 ("7.38% deficit flexible vs rigid") uses the Frontiers/Zhou/Ma operating
point (V=10.59 m/s, Ω=7.55 rpm), while V-01 validation uses NREL rated
(V=10.659 m/s, Ω=7.518 rpm). The two points differ; without cross-evaluation
the 7.38% claim could carry contamination from the BEM bias (+3.34% T,
+2.67% Q already known from V-01).

This script evaluates the SAME BEMSolver+blade at BOTH operating points and
reports the deltas so the aeroelastic deficit is not mixed with operating-point
differences.

Output: docs/validation_data/generated/cross_baseline_rigid_bem.csv
        docs/validation_data/generated/cross_baseline_rigid_summary.md
"""

from __future__ import annotations

import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
BLADE_YAML = REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"
OUT_DIR = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Operating points to evaluate.
POINTS = [
    {"label": "V-01 (NREL rated)", "v_inf": 10.659, "omega_rpm": 7.518, "pitch_deg": 0.0},
    {"label": "V-04 (Frontiers)", "v_inf": 10.59, "omega_rpm": 7.55, "pitch_deg": 0.0},
]

COMMON = {
    "yaw_deg": 0.0,
    "rho": 1.225,
    "mu": 1.81206e-5,
    "hub_height": 150.0,
    # The FSI campaign yamls (frontiersin_2025_yaw) use shear_exp: 0.0 —
    # a shear-based rigid anchor misaligned with the campaign inflates the
    # flexible-vs-rigid comparison by ~1.7% (see Engram: both_off closure).
    "shear_exp": 0.0,
}

# Reference values reported in the dossier (for sanity check).
DOSSIER_RIG_TORQUE = 20.081e6  # N·m at the V-04 point per plot_torque_signal_analysis.py
DOSSIER_FLEX_TORQUE = 18.599e6  # N·m corotational FSI at V-04 point


def main() -> int:
    if not BLADE_YAML.exists():
        print(f"ERROR: blade YAML not found at {BLADE_YAML}", file=sys.stderr)
        return 2

    from aeroelast.models.blade.aerodynamics import load_blade_aero
    from aeroelast.solvers.bem.engine import BEMSolver

    print(f"[A7] Loading blade aero from {BLADE_YAML.name}...")
    blade = load_blade_aero(str(BLADE_YAML))

    results = []
    for pt in POINTS:
        bem = BEMSolver(
            blade,
            rho=COMMON["rho"],
            mu=COMMON["mu"],
            yaw=COMMON["yaw_deg"],
            hub_height=COMMON["hub_height"],
            shear_exp=COMMON["shear_exp"],
        )
        print(f"[A7] Evaluating {pt['label']}: V={pt['v_inf']} m/s, Ω={pt['omega_rpm']} RPM...")
        out = bem.compute(v_inf=pt["v_inf"], omega=pt["omega_rpm"], pitch=pt["pitch_deg"])
        results.append({
            "point": pt["label"],
            "v_inf_m_s": pt["v_inf"],
            "omega_rpm": pt["omega_rpm"],
            "pitch_deg": pt["pitch_deg"],
            "torque_Nm": float(out.torque),
            "power_W": float(out.power),
            "thrust_N": float(out.thrust),
            "C_T": float(out.CT) if hasattr(out, "CT") else float("nan"),
            "C_P": float(out.CP) if hasattr(out, "CP") else float("nan"),
        })

    # Export CSV.
    csv_path = OUT_DIR / "cross_baseline_rigid_bem.csv"
    keys = list(results[0].keys())
    with csv_path.open("w") as f:
        f.write(",".join(keys) + "\n")
        for r in results:
            f.write(
                ",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys)
                + "\n"
            )
    print(f"[A7] Wrote {csv_path}")

    # Summary markdown.
    md = [
        "# A7 — Cross-baseline BEM-only at V-01 and V-04 operating points",
        "",
        "Same `BEMSolver` + same blade YAML (`IEA-15-240-RWT.yaml`) evaluated at both",
        "the V-01 (NREL rated) and V-04 (Frontiers) operating points. This isolates",
        "the operating-point contribution from the flexible-vs-rigid deficit.",
        "",
        "| Operating point | V_inf [m/s] | Ω [RPM] | Torque [MN·m] | Power [MW] | Thrust [MN] |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        md.append(
            f"| {r['point']} | {r['v_inf_m_s']:.3f} | {r['omega_rpm']:.3f} | "
            f"{r['torque_Nm'] / 1e6:.3f} | {r['power_W'] / 1e6:.3f} | {r['thrust_N'] / 1e6:.3f} |"
        )

    # Comparisons.
    q_v01 = results[0]["torque_Nm"]
    q_v04 = results[1]["torque_Nm"]
    delta_q = (q_v04 - q_v01) / q_v01 * 100.0
    md += [
        "",
        f"**Delta torque V-04 vs V-01 (rigid BEM only)**: {delta_q:+.2f}%",
        f"({(q_v04 - q_v01) / 1e6:+.3f} MN·m).",
        "",
        '## Decomposition of the §3.7.1 "flexible vs rigid" claim',
        "",
        f"- Dossier rigid baseline (used in §3.7.1, V-04 point): {DOSSIER_RIG_TORQUE / 1e6:.3f} MN·m",
        f"- This run rigid at V-04 point: {q_v04 / 1e6:.3f} MN·m",
        f"- This run rigid at V-01 point: {q_v01 / 1e6:.3f} MN·m",
        f"- Dossier FSI flexible (corotational at V-04 point): {DOSSIER_FLEX_TORQUE / 1e6:.3f} MN·m",
        "",
        f"Deficit flexible vs rigid (V-04 point, §3.7.1 number): "
        f"{(DOSSIER_FLEX_TORQUE - DOSSIER_RIG_TORQUE) / DOSSIER_RIG_TORQUE * 100:+.2f}%",
        f"Deficit flexible vs rigid (V-01 point, hypothetical): "
        f"{(DOSSIER_FLEX_TORQUE - q_v01) / q_v01 * 100:+.2f}%",
        "",
        f"Operating-point contribution relative to the V-01 rigid baseline: {delta_q:+.2f}%.",
        "This is why the aeroelastic deficit must be computed against the V-04 rigid",
        "baseline, not against the NREL-rated V-01 rigid point.",
    ]
    md_path = OUT_DIR / "cross_baseline_rigid_summary.md"
    md_path.write_text("\n".join(md) + "\n")
    print(f"[A7] Wrote {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
