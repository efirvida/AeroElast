# Chapter 6 DLC Prepared Cases

Prepared case templates for Chapter 6 DLC runs.
The solver YAMLs are now solver-only (no dlc metadata block).
Use ch6_run_matrix.csv as the execution plan source.

Reference simulation bases:
- /scratch/leahk/eduardo.donestevez/simulations/bem_0_10 (corotational + fluid)
- /scratch/leahk/eduardo.donestevez/simulations/bem_0_10_full (inertial when available)

Files:
- fluid_dlc_*.yaml
- solid_corotational_dlc_*.yaml
- solid_inertial_dlc_*.yaml
- ch6_run_matrix.csv (full execution matrix)

Expected simulation counts per solver (from Table 6-1):
- DLC 1.1: 72
- DLC 1.3: 72
- DLC 1.4: 6
- DLC 1.5: 48
- DLC 6.1: 12
- DLC 6.3: 12
- Total per solver: 222
- Total both solvers: 444

How to use:
- Pick template pair (solid_* + fluid_*) for desired DLC and solver
- Iterate the cases listed in ch6_run_matrix.csv
- For each run, inject wind/event/yaw/seed from the matrix into your launcher/runtime args
