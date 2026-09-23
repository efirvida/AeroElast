"""Campbell diagram: rotating modal frequencies vs rotor speed (shell).

Sweeps Omega over [0, rated+] and records the 1st flap / 1st edge / 2nd flap
frequencies from the rotating-frame modal (K + K_G + K_SP), then plots the
diagram with the 1P/2P/3P/4P per-rev lines and the reference points.

Usage:
  python tools/run_sx_campbell.py --rpms 0 3 5 7.56 9 11 \
      --csv docs/validation_data/generated/sx_campbell.csv \
      --plot docs/validation_plots/figures/campbell.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))


def shell_campbell(rpms, element_size=0.5):
    from run_s4_rotating_modal import assemble_rotating, classify_modes

    asm, mesh = None, None
    rows = []
    for rpm in rpms:
        if asm is None:
            from run_s5_oneway import build_assembly
            asm, mesh = build_assembly(element_size)
        omega = rpm * 2 * np.pi / 60.0
        freqs, modes, free = assemble_rotating(asm, mesh, omega, n_modes=10)
        labels = classify_modes(mesh, modes, free)
        flap = sorted(f for f, l in zip(freqs, labels) if l == "flap")
        edge = sorted(f for f, l in zip(freqs, labels) if l == "edge")
        rows.append({
            "rpm": rpm,
            "f1f": flap[0] if flap else np.nan,
            "f2f": flap[1] if len(flap) > 1 else np.nan,
            "f1e": edge[0] if edge else np.nan,
        })
        print(f"rpm={rpm:5.1f}: 1F={rows[-1]['f1f']:.4f}  2F={rows[-1]['f2f']:.4f}  "
              f"1E={rows[-1]['f1e']:.4f} Hz", flush=True)
    return rows


def plot_campbell(rows, out_path, of_ref=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rpm = np.array([r["rpm"] for r in rows])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(rpm, [r["f1f"] for r in rows], "o-", label="1st flap (shell)")
    ax.plot(rpm, [r["f2f"] for r in rows], "s-", label="2nd flap (shell)")
    ax.plot(rpm, [r["f1e"] for r in rows], "^-", label="1st edge (shell)")

    for p in (1, 2, 3, 4):
        ax.plot(rpm, rpm * p / 60.0, "--", color="gray", lw=0.7,
                label=f"{p}P" if p == 1 else None)

    if of_ref is not None:
        pts = of_ref.get("all", [])
        if pts:
            r, f = zip(*pts)
            ax.plot(r, f, "x", color="tab:orange", ms=6, mfc="none",
                    label="OpenFAST MBC3")

    ax.set_xlabel("rotor speed (rpm)")
    ax.set_ylabel("frequency (Hz)")
    ax.set_title("Campbell diagram — IEA 15 MW blade (rotating-frame modal)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"plot written to {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpms", type=float, nargs="+", default=[0, 3, 5, 7.56, 9, 11])
    ap.add_argument("--csv", type=Path,
                    default=Path("docs/validation_data/generated/sx_campbell.csv"))
    ap.add_argument("--plot", type=Path,
                    default=Path("docs/validation_plots/figures/campbell.png"))
    ap.add_argument("--element-size", type=float, default=0.5)
    args = ap.parse_args()

    rows = shell_campbell(args.rpms, args.element_size)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with open(args.csv, "w") as fh:
            fh.write("rpm,f1f_Hz,f2f_Hz,f1e_Hz\n")
            for r in rows:
                fh.write(f"{r['rpm']},{r['f1f']:.6f},{r['f2f']:.6f},{r['f1e']:.6f}\n")
        print(f"wrote {args.csv}")

    # overlay the OpenFAST MBC3 reference points when available
    of_ref = None
    of_csv = Path("docs/validation_data/generated/sx_campbell_openfast.csv")
    if of_csv.exists():
        d = np.loadtxt(of_csv, delimiter=",", skiprows=1)
        if d.ndim == 1:
            d = d[None, :]
        of_ref = {"all": [(r, f) for r, f in zip(d[:, 0], d[:, 1])]}
    plot_campbell(rows, args.plot, of_ref=of_ref)


if __name__ == "__main__":
    main()
