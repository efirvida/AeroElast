# Aero-FSI Report Contract

Generic Aero-FSI runs now emit a small schema manifest plus up to three CSV reports in the output folder. This gives post-processing and validation scripts a stable contract without having to inspect the participant code.

## Quick Path

1. Open `aero_report_schema.json` in the output folder.
2. Pick the report file that matches the granularity you need.
3. Use the `fieldnames` array from the manifest instead of inferring column meaning from the first CSV row.

From the terminal you can inspect the manifest quickly with:

```bash
fem-shell-aero-report-info path/to/output
```

If you want one report only, or its ordered fields, use:

```bash
fem-shell-aero-report-info path/to/output --report aero_sectional_report
fem-shell-aero-report-info path/to/output --report aero_sectional_report --fields
```

## Report Map

| File | When it appears | What it represents |
|------|------------------|--------------------|
| `aero_report_schema.json` | whenever any Aero-FSI report is written | Stable manifest with schema versions, field order, and report summaries |
| `aero_report.csv` | every converged Aero-FSI time window | Rotor-global loads plus blade-integrated diagnostics |
| `aero_spanwise_report.csv` | when backend metadata exposes per-panel data | One row per aerodynamic panel |
| `aero_sectional_report.csv` | when backend metadata exposes per-panel data | One row per blade radial bin aggregated from panels |

## Start Here

| If you need | Read |
|-------------|------|
| Total rotor load history | `aero_report.csv` |
| Blade-by-blade integrated forces and moments | `aero_report.csv` |
| Per-panel local wind and panel force | `aero_spanwise_report.csv` |
| Validation-friendly radial trends without panel-level noise | `aero_sectional_report.csv` |
| Stable field order for parsers | `aero_report_schema.json` |

## Key Coordinates

| Column | Meaning |
|--------|---------|
| `Radius [m]` | Absolute panel or bin radius measured from the rotor center |
| `Radius Fraction [-]` | Panel radius normalized by the blade tip radius |
| `Spanwise [m]` | Absolute spanwise coordinate within the blade aerodynamic surface |
| `Blade Span Fraction [-]` | Panel spanwise coordinate normalized within that blade |
| `Radius Center Fraction [-]` | Sectional-bin center radius normalized by blade tip radius |
| `Blade Span Center Fraction [-]` | Sectional-bin center span normalized within that blade |

## Manifest Shape

The schema manifest contains these top-level keys:

| Key | Meaning |
|-----|---------|
| `manifest_version` | Version of the JSON manifest format itself |
| `backend` | Active aerodynamic backend, for example `bem` or `vlm` |
| `n_blades` | Blade count represented in the report files |
| `sectional_bins` | Configured radial-bin count for the sectional report |
| `reports` | Per-report metadata keyed by report type |

Each report entry includes:

| Key | Meaning |
|-----|---------|
| `file` | CSV filename in the output folder |
| `schema_version` | Version of that CSV schema |
| `granularity` | Row meaning, such as per-window, per-panel, or per-bin |
| `summary` | Short human-readable description |
| `fieldnames` | Stable ordered list of CSV columns |

## Checklist

- [ ] My parser reads `aero_report_schema.json` before consuming any CSV.
- [ ] My Python loader uses `load_aero_report_table(...)` instead of hard-coded CSV assumptions.
- [ ] I use `schema_version` and `fieldnames` from the manifest, not ad hoc header assumptions.
- [ ] I choose `aero_report.csv` for integrated history, `aero_spanwise_report.csv` for panel detail, or `aero_sectional_report.csv` for radial trends.

## Next Step

If you need plots or validators, build them against `aero_report_schema.json` first and then map the specific CSV fields you care about. That keeps downstream tooling stable even if optional columns evolve later.