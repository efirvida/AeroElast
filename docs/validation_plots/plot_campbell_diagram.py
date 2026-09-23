"""A1 — Campbell diagram (static level 1).

Plots the IEA 15 MW blade natural frequencies (from V-02, computed offline by
``compute_modal_participation.py``) against the harmonic excitation lines mP
(m = 1, 3, 6, 9) over the rotor speed range Ω ∈ [0, 12] RPM. Highlights the
operating point Ω_rated = 7.518 RPM and any modal-harmonic crossing inside the
operating envelope.

LEVEL 1 (this script): static eigenvalues — does NOT include centrifugal
stiffening K_G(Ω) or spin softening K_SP(Ω). Valid as a first cut for
diagnosing resonance proximity; the rotational correction shifts flapwise
frequencies up by ~3-8% near rated for blades of this scale.

LEVEL 2 (TODO): the ModalSolver currently reads K/M from internal Rust hooks
and does not accept an externally-assembled K_eff = K + K_G(Ω) + K_SP(Ω). To
produce a rotational Campbell diagram, add an overload to
``modal_solve_coo`` accepting external (rows, cols, vals) triplets for K, and
loop over Ω rebuilding K_G via ``domain.assemble_geometric_stiffness(omega)``
and K_SP analytically. Out of scope for this script.

Reads:
  docs/validation_data/generated/v02_modal_participation.csv (from A2)

Writes:
  docs/validation_plots/figures/fig_v02_campbell_diagram.png
  docs/validation_plots/figures/fig_v02_campbell_diagram.pdf
  docs/validation_data/generated/campbell_crossings.csv
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA_IN = REPO_ROOT / "docs" / "validation_data" / "generated" / "v02_modal_participation.csv"
OUT_FIG = REPO_ROOT / "docs" / "validation_plots" / "figures"
OUT_DATA = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_DATA.mkdir(parents=True, exist_ok=True)

OMEGA_RATED_RPM = 7.518
OMEGA_MAX_RPM = 12.0
HARMONICS = [1, 3, 6, 9]   # 3 blades → 3P, 6P, 9P are the dominant excitations.
MAX_FREQ_HZ = 5.0          # cap y-axis to keep the diagram readable.


def label_mode(row: pd.Series) -> str:
    """Identify the dominant character of a mode from the participation row.

    Convention (span = global Z, root clamped):
      Rx-dominant → flapwise bending
      Ry-dominant → edgewise bending
      Rz-dominant → torsion
      Tx-dominant → 2nd+ edgewise (where translation dominates over rotation)
      Ty-dominant → 2nd+ flapwise
    """
    fields = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]
    values = {f: row[f] for f in fields}
    dom = max(values, key=values.get)
    if dom in ("Rx", "Ty"):
        return "flapwise"
    if dom in ("Ry", "Tx"):
        return "edgewise"
    if dom == "Rz":
        return "torsion"
    return "axial"


def main() -> int:
    if not DATA_IN.exists():
        print(f"ERROR: missing {DATA_IN} — run compute_modal_participation.py first.")
        return 2

    df = pd.read_csv(DATA_IN)
    df["character"] = df.apply(label_mode, axis=1)
    print(df[["mode", "freq_hz", "character"]].to_string(index=False))

    # Group modes by character to assign labels and colors.
    char_colors = {"flapwise": "#c0392b", "edgewise": "#2980b9", "torsion": "#27ae60", "axial": "#7f8c8d"}
    char_counter = {k: 0 for k in char_colors}
    mode_labels = []
    for _, row in df.iterrows():
        c = row["character"]
        char_counter[c] += 1
        n = char_counter[c]
        prefix = {"flapwise": "F", "edgewise": "E", "torsion": "T", "axial": "A"}[c]
        mode_labels.append(f"{n}{prefix}  ({row['freq_hz']:.3f} Hz)")

    # ── Plot ──────────────────────────────────────────────────────────────
    omega = np.linspace(0, OMEGA_MAX_RPM, 200)
    omega_hz = omega / 60.0

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    # Harmonic excitation lines mP.
    harmonic_color = "#34495e"
    for m in HARMONICS:
        ax.plot(omega, m * omega_hz, color=harmonic_color, ls="--", lw=0.9, alpha=0.7)
        # Label at the right edge of the plot.
        y_end = m * omega_hz[-1]
        if y_end < MAX_FREQ_HZ:
            ax.text(omega[-1] * 1.005, y_end, f"{m}P", color=harmonic_color, fontsize=8,
                    va="center", ha="left")

    # Horizontal lines for natural frequencies (static — no Ω dependence yet).
    crossings = []
    for i, row in df.iterrows():
        if row["freq_hz"] > MAX_FREQ_HZ:
            continue
        color = char_colors[row["character"]]
        ax.axhline(row["freq_hz"], color=color, lw=1.2, alpha=0.85, zorder=2)
        ax.text(0.15, row["freq_hz"] + 0.04, mode_labels[i], color=color, fontsize=8)
        # Detect crossings with mP lines inside [0, OMEGA_MAX_RPM].
        for m in HARMONICS:
            # mP(Ω) = m·Ω/60. Crossing when m·Ω/60 = f_n → Ω = 60·f_n/m.
            omega_cross = 60.0 * row["freq_hz"] / m
            if 0 < omega_cross <= OMEGA_MAX_RPM:
                margin_to_rated = (omega_cross - OMEGA_RATED_RPM) / OMEGA_RATED_RPM
                crossings.append({
                    "mode": int(row["mode"]),
                    "character": row["character"],
                    "freq_hz": float(row["freq_hz"]),
                    "harmonic": m,
                    "omega_crossing_rpm": float(omega_cross),
                    "margin_to_rated_pct": float(100.0 * margin_to_rated),
                })

    # Operating point marker.
    ax.axvline(OMEGA_RATED_RPM, color="black", lw=1.5, ls=":", alpha=0.9)
    ax.text(OMEGA_RATED_RPM + 0.1, MAX_FREQ_HZ * 0.97,
            f"Ω_rated = {OMEGA_RATED_RPM:.3f} RPM",
            color="black", fontsize=9, va="top", rotation=0)

    ax.set_xlim(0, OMEGA_MAX_RPM)
    ax.set_ylim(0, MAX_FREQ_HZ)
    ax.set_xlabel("Rotor speed Ω [RPM]")
    ax.set_ylabel("Frequency [Hz]")
    ax.set_title("Campbell diagram — IEA 15 MW blade (static, no centrifugal stiffening)\n"
                 "Natural frequencies (horizontal) vs. harmonic excitation mP (dashed)")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()

    png = OUT_FIG / "fig_v02_campbell_diagram.png"
    pdf = OUT_FIG / "fig_v02_campbell_diagram.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(f"[A1] Saved {png}")

    # Export crossing table.
    if crossings:
        cdf = pd.DataFrame(crossings).sort_values("omega_crossing_rpm")
        cdf.to_csv(OUT_DATA / "campbell_crossings.csv", index=False)
        print(f"[A1] Wrote {OUT_DATA / 'campbell_crossings.csv'} ({len(cdf)} crossings)")

        # Identify any crossing within ±5% of rated.
        near_rated = cdf[abs(cdf["margin_to_rated_pct"]) < 5.0]
        if len(near_rated) > 0:
            print("[A1] WARNING — modal/harmonic crossings within ±5% of rated:")
            print(near_rated.to_string(index=False))
        else:
            print(f"[A1] No modal/harmonic crossings within ±5% of Ω_rated={OMEGA_RATED_RPM} RPM.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
