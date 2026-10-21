"""B4 K1/K4/K5 measurement instrument (no verdicts or tolerance changes).

Run with --dt-scale {1,0.5} --output PATH; stdout is the plain measurement table.
JSON stores fractions (not percentages), angles in degrees, and exact control prints.
Finding: write_aerodyn_dvr formats dt with two decimals. Scale 1 uses its output
byte-for-byte; scale 0.5 rewrites only combined-case dt tokens to preserve 0.125.
"""

import argparse
import contextlib
import io
import json
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.validation.bem import test_bem_openfast_parity as p  # noqa: E402


def solver(blade, yaw=0.0, shear=0.0):
    return p.BEMSolver(
        blade, rho=1.225, mu=1.464e-5, precone=p.ob.IEA15MW_PRECONE_DEG,
        tilt=p.ob.IEA15MW_SHAFT_TILT_DEG, yaw=yaw,
        hub_height=p.ob.IEA15MW_HUB_HEIGHT, shear_exp=shear,
    )


def reference(cases, dt, tmax, scale):
    driver = p._openfast_driver_or_skip()
    with tempfile.TemporaryDirectory(prefix="b4-t7-") as td:
        td = Path(td)
        primary = p.ob.write_bem_primary(p.PRIMARY, td / "AD_clean.dat")
        blade = p.ob.build_blade_aero_from_aerodyn(primary)
        dvr = p.ob.write_aerodyn_dvr(td / "measure.dvr", primary, cases, dt=dt, tmax=tmax)
        if scale == 0.5:
            lines = dvr.read_text().splitlines(keepends=True)
            start = next(i for i, line in enumerate(lines) if line.startswith("HWndSpeed")) + 2
            for i in range(start, start + len(cases)):
                tokens = lines[i].split()
                tokens[5] = str(dt)
                lines[i] = " ".join(tokens) + "\n"
            dvr.write_text("".join(lines))
        refs = p.ob.parse_aerodyn_outs(p.ob.run_aerodyn_driver(driver, dvr, td))
    assert len(refs) == len(cases)
    return blade, refs


def measure(scale):
    blade, refs = reference(p.CASES, scale, 60.0, scale)
    yaw_blade, yaw_refs = reference(p.A2_CASES, p.A2_DT * scale, p.A2_TMAX, scale)
    r = p.ob.IEA15MW_HUB_RAD + p.ob.parse_aerodyn_blade(p.BLADE).bl_spn
    mask = (r > p.INTERIOR[0] * r[-1]) & (r < p.INTERIOR[1] * r[-1])
    parity = dict(r=r, solver=solver(blade), reference=refs, interior=mask)
    yaw_data = dict(blade_aero=yaw_blade, reference=yaw_refs)
    yaml = p.load_blade_aero(str(p.YAML_BLADE))
    rows, control = {}, io.StringIO()
    for label, cases, references, aero in (
        ("A1", p.CASES, refs, blade), ("A2", p.A2_CASES, yaw_refs, yaw_blade),
        ("A3", p.CASES, refs, yaml),
    ):
        for i, (case, ref) in enumerate(zip(cases, references)):
            wind, rpm, pitch = case[:3]
            shear, yaw = case[3:] if label == "A2" else (0.0, 0.0)
            computed = solver(aero, yaw, shear).compute(wind, rpm, pitch)
            rows[f"{label}.{i}"] = dict(
                thrust=(computed.thrust - ref.thrust) / ref.thrust,
                torque=(computed.torque - ref.torque) / ref.torque,
            )
            with contextlib.redirect_stdout(control):
                if label == "A1":
                    p.test_rotor_performance_matches_aerodyn(parity, i)
                    p.test_spanwise_loads_match_aerodyn(parity, i)
                elif label == "A2":
                    p.test_yaw_and_shear_match_aerodyn(yaw_data, i)
                else:
                    p.test_bem_with_repo_default_polars_matches_aerodyn(parity, yaml, i)
            if label == "A1":
                da = float(np.abs(computed.alpha - ref.alpha_deg)[mask].max())
                ar = float(np.abs(ref.alpha_deg)[mask].max())
                cn = np.abs((computed.cn - ref.cn) / np.where(ref.cn == 0, np.nan, ref.cn))
                rows[f"spanwise.{i}"] = dict(
                    max_dalpha_deg=da, max_ref_alpha_deg=ar,
                    e_alpha=da / ar, mean_abs_dcn=float(np.nanmean(cn[mask])),
                )
    viterna = []
    with contextlib.redirect_stdout(control), patch.object(
        p, "assert_residual_below", side_effect=lambda value, **kw: viterna.append(value)
    ):
        p.test_viterna_extension_matches_official_from_official_data()
    rows["Viterna"] = dict(dCl=viterna[0], dCd=viterna[1])
    return dict(dt_scale=scale, a1_dt=scale, a2_dt=p.A2_DT * scale,
                rows=rows, control_lines=control.getvalue().splitlines())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dt-scale", type=float, choices=(1.0, 0.5), required=True)
    parser.add_argument("--output", type=Path, required=True, help="JSON output file")
    parser.add_argument("--control", type=Path, help="require every printed metric in T6 log")
    args = parser.parse_args()
    data = measure(args.dt_scale)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    print(f"dt-scale={args.dt_scale}; A1 dt={data['a1_dt']}; A2 dt={data['a2_dt']}")
    print("row | metric=value (fractions unless _deg; full float representation)")
    for row, metrics in data["rows"].items():
        print(row + " | " + " | ".join(f"{k}={v!r}" for k, v in metrics.items()))
    if args.control:
        baseline = args.control.read_text()
        missing = [line for line in data["control_lines"] if line not in baseline]
        print(f"T6 control: {len(data['control_lines']) - len(missing)}/"
              f"{len(data['control_lines'])} printed lines reproduced")
        if missing:
            raise SystemExit("T6 CONTROL MISMATCH:\n" + "\n".join(missing))


if __name__ == "__main__":
    main()
