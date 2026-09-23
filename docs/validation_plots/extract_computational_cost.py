"""
Computational cost extraction from solver profiling logs.

The AeroElast Rust solver emits per-window profiling lines of the form

    RotorFsi window N t=...s wall=...ms, checkpoint_save: a (×k_1),
    precice_read: b (×k_2), ... kg_update_if_needed: x (×k_n)

This script parses ``solid.log`` from one or more campaigns and produces a
summary table of wall-clock time per phase, normalised per simulated second.
The result feeds §5.6 (computational cost) of the validation document.

Outputs
-------
docs/figures/fig_5_6_cost_breakdown.png  (300 dpi)
docs/figures/fig_5_6_cost_breakdown.pdf
stdout: a markdown-formatted summary table ready to paste into the doc.
"""
import pathlib
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
LOGS = [
    ("bem_0_10",     BASE / "simulations" / "bem_0_10"     / "solid.log"),
    ("bem_90_50_S",  BASE / "simulations" / "bem_90_50_S"  / "solid.log"),
]
OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Regex: capture the window number, total wall time, and each phase entry.
LINE_RE = re.compile(
    r"RotorFsi window (\d+) t=([\d.]+)s wall=([\d.]+)ms, (.+)"
)
PHASE_RE = re.compile(r"(\w+): ([\d.]+)ms \(×(\d+)\)")


def parse_log(path: pathlib.Path) -> dict[str, dict[str, list[float]]]:
    """Return dict {phase: [ms_per_window, ...]}, plus 'wall' total per window."""
    stats = defaultdict(list)
    if not path.exists():
        return {}
    with open(path, errors="ignore") as f:
        for line in f:
            m = LINE_RE.search(line)
            if not m:
                continue
            window = int(m.group(1))
            sim_time = float(m.group(2))
            wall_ms = float(m.group(3))
            stats["__wall_total__"].append(wall_ms)
            stats["__window__"].append(window)
            stats["__sim_time__"].append(sim_time)
            for pm in PHASE_RE.finditer(m.group(4)):
                name = pm.group(1)
                ms = float(pm.group(2))
                count = int(pm.group(3))
                # Total wall ms spent on this phase in this window
                stats[name].append(ms * count)
    return stats


# ---------------------------------------------------------------------------
# Parse all logs
# ---------------------------------------------------------------------------
summary = {}
for tag, path in LOGS:
    print(f"\nParsing {tag}: {path}")
    stats = parse_log(path)
    if not stats:
        print("  (no profiling lines found — skipping)")
        continue
    nw = len(stats["__wall_total__"])
    print(f"  {nw} windows profiled")
    summary[tag] = stats

if not summary:
    raise RuntimeError("No profiling data found in any log")

# ---------------------------------------------------------------------------
# Aggregate per-phase totals and print Markdown table
# ---------------------------------------------------------------------------
print("\n" + "=" * 90)
print("Computational cost breakdown — wall time per simulated second")
print("=" * 90)
header = "| Run | Δt [s] | N_win | wall/sim [s/s] | Newmark | K_G upd | preCICE | callback | other |"
print(header)
print("|" + "|".join("---" for _ in header.split("|")[1:-1]) + "|")

for tag, stats in summary.items():
    n_win = len(stats["__wall_total__"])
    total_wall_s = sum(stats["__wall_total__"]) / 1000.0
    sim_times = sorted(stats["__sim_time__"])
    dt = sim_times[1] - sim_times[0] if len(sim_times) > 1 else 0.0
    sim_total = sim_times[-1] - sim_times[0] + dt
    wall_per_sim = total_wall_s / sim_total if sim_total > 0 else 0.0

    def frac(phase: str) -> float:
        if phase not in stats:
            return 0.0
        return 100.0 * sum(stats[phase]) / sum(stats["__wall_total__"])

    f_newmark = frac("newmark_step")
    f_kg = frac("kg_update_if_needed")
    f_precice = frac("precice_advance")
    f_callback = frac("callback")
    f_other = 100 - (f_newmark + f_kg + f_precice + f_callback)

    print(f"| {tag} | {dt:.4f} | {n_win} | {wall_per_sim:.2f} | "
          f"{f_newmark:.1f} % | {f_kg:.1f} % | {f_precice:.1f} % | "
          f"{f_callback:.1f} % | {f_other:.1f} % |")

print()

# ---------------------------------------------------------------------------
# Stacked-bar plot per run
# ---------------------------------------------------------------------------
runs = list(summary.keys())
phases = ["newmark_step", "kg_update_if_needed", "precice_advance",
          "callback", "checkpoint_save", "checkpoint_restore",
          "force_transform", "inertial_forces", "disp_gather_transform",
          "precice_read", "precice_write", "force_preproc",
          "kt_update_if_needed"]
phase_colors = plt.get_cmap("tab20")(np.linspace(0, 1, len(phases)))

fractions = np.zeros((len(phases), len(runs)))
for j, tag in enumerate(runs):
    stats = summary[tag]
    total = sum(stats["__wall_total__"])
    for i, ph in enumerate(phases):
        if ph in stats:
            fractions[i, j] = 100.0 * sum(stats[ph]) / total

fig, ax = plt.subplots(figsize=(9, 5.5))
bottom = np.zeros(len(runs))
for i, ph in enumerate(phases):
    if fractions[i].max() < 0.1:
        continue   # skip phases under 0.1 % everywhere
    ax.bar(runs, fractions[i], bottom=bottom,
           color=phase_colors[i], label=ph)
    bottom += fractions[i]

ax.set_ylabel("Wall-time share [%]")
ax.set_title("AeroElast computational cost breakdown by phase\n"
             "(per-window profiling, wall time × invocations / window)")
ax.set_ylim(0, 100)
ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", fontsize=8, frameon=False)
ax.grid(axis="y", alpha=0.3)

fig.tight_layout()
png = OUT_DIR / "fig_5_6_cost_breakdown.png"
pdf = OUT_DIR / "fig_5_6_cost_breakdown.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
