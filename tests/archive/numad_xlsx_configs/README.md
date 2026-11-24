# numad_xlsx_configs

Case configs that built the IEA-15 blade from the UTD NuMAD spreadsheet (`excel_file`), archived by
gmsh-official-mesher T6b (S6, user decision L-S41 of `odd/tasks/gmsh-official-mesher.md`).

The NuMAD mesher is retired (S5) and `BladeMesh` / `RotorMesh` reject `excel_file` (S7a): these configs no longer load.  The
UTD blade has a thinner outboard trailing-edge reinforcement than the official windIO input (AGENTS.md rule 7,
`docs/blade_input_divergence_utd_vs_official.md`), so none of them is a current anchor.  They are kept as a record of the
input-sensitivity runs only; their relative paths (`../../...`) are relative to their former directory under `tests/`.
