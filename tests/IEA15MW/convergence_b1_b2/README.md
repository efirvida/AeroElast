# B1 + B2 — h/dt convergence campaign for V-04 corotacional FSI

**Goal**: close F5 (the only remaining blocking gap of the integral review for paper P1).
Without these runs, the manuscript cannot pass first-round review at any journal of
computational mechanics.

## Scope

Single operating point — yaw=0°, corotacional, V=10.59 m/s, Ω=7.55 RPM, pitch=0°.
This is the same nominal point as `frontiersin_2025_yaw/solid_corotational_yaw_0.yaml`,
which serves as the baseline. The baseline is **NOT re-run** — its outputs already exist
in `$SCRATCH/frontiersin_results_corotational/yaw_0/`.

Five points total, four new runs:

| Tag | element_size [m] | time-window-size [s] | n_samples | Notes |
|---|---:|---:|---:|---|
| `h_coarse`     | 0.500 | 0.01 | 300 | B1 coarse (×2) |
| `h_baseline`   | 0.250 | 0.01 | 300 | already in frontiersin baseline |
| `h_fine`       | 0.125 | 0.01 | 300 | B1 fine (÷2) |
| `dt_coarse`    | 0.250 | 0.02 | 300 | B2 coarse (×2) |
| `dt_baseline`  | 0.250 | 0.01 | 300 | same baseline |
| `dt_fine`      | 0.250 | 0.005 | 300 | B2 fine (÷2) |

Simulated time: **30 s** per case (sufficient to reach steady state past 20 s window).

Each run uses its dedicated `precice-config-*.xml` (only `time-window-size` differs)
and its own `solid_*.yaml` / `fluid_*.yaml`.

**B1 design (corrected 2026-09-08)**: the fluid mesh is PINNED at `element_size =
0.25` for all B1 levels, so the BEM station radii are identical across h levels
and B1 isolates structural h-convergence only. The first B1 runs (2026-05-23)
co-varied the fluid `element_size` with the solid mesh, which shifted the BEM
stations (r0 = 0.468 / 0.554 / 0.585 m) and produced a non-monotonic flap_mean
(negative Richardson order) — the confounded campaign directory was removed
to free quota space; the clean re-run (`runB1_rerun.srm`) writes to
`$SCRATCH/convergence_b1_b2_results_b1fix/`.

**B2 re-run (2026-09-08)**: the original `dt_coarse`/`dt_fine` results were
removed together with the confounded campaign directory, and the raw
`bem_report.csv` files are required to regenerate the convergence figure.
`runB1_rerun.srm` therefore re-runs all four cases (h + dt).

**Disk footprint (2026-09-08)**: the first re-run attempt died with
`OSError: [Errno 122] Disk quota exceeded` (group `leahk` hit its 2 TiB hard
limit) because the solid checkpoints wrote `fields.vtu` + `deformed_mesh.h5`
every 0.05 s (~94 MB/checkpoint on h_fine ≈ 56 GB per 30 s run). All four solid
YAMLs now use `write_interval: 0.5` + `save_deformed_mesh: false`: the B1/B2
metrics come exclusively from `fluid/bem_report.csv`, and `state.npz`
checkpoints (restart capability, ≤0.5 s lost) are kept. Total campaign
footprint ≈ 9 GB instead of ≈ 100 GB.

## Outputs

By default to `$SCRATCH/convergence_b1_b2_results_b1fix/<tag>/`:

```
<tag>/
  corotational/   # state.npz checkpoints + fields.vtu every 0.5 s (reduced footprint)
  fluid/          # bem_report.csv + bem_sectional.csv per timestep
  logs/
    solid.log
    fluid.log
```

## Postprocessing (after all runs complete)

`docs/validation_plots/plot_convergence_b1_b2.py` plots:
- Flap mean / std / p-p vs element_size (B1)
- Flap mean / std / p-p vs time-window-size (B2)
- Power mean / std vs both
- Richardson extrapolation estimate of the asymptotic value (order of convergence ~2)

```bash
RESULTS_BASE=$SCRATCH/convergence_b1_b2_results_b1fix \
BASELINE_PATH=$SCRATCH/frontiersin_results_corotational_100s/yaw_0 \
python docs/validation_plots/plot_convergence_b1_b2.py
```

Outputs: `docs/validation_data/generated/convergence_b1_b2_summary.{csv,md}` and
`docs/validation_plots/figures/fig_v04_convergence_b1_b2.{png,pdf}`.

## Submission

```bash
cd tests/IEA15MW/convergence_b1_b2
sbatch runB1_rerun.srm
```

Override outputs base via env var:
```bash
RESULTS_BASE=$SCRATCH/my_runs sbatch runB1_rerun.srm
```

## Acceptance criteria (per dossier §1 grietas)

The campaign closes F5 when:
1. B1: |flap_mean(h_fine) − flap_mean(h_baseline)| / flap_mean(h_baseline) < 1 %
   AND |flap_mean(h_baseline) − flap_mean(h_coarse)| / flap_mean(h_baseline) < 5 %
   (i.e. baseline is asymptotic enough)
2. B2: same criterion replacing element_size with time-window-size
3. Both: solver completes without divergence, residuals < 1×10⁻⁴ as in §5.6

If a finer level diverges or shows non-monotonic behaviour, that is itself a result —
it bounds the usable parameter range and must be reported in the §V-04 robustness
sub-section.
