# T31 finding: no fibre-orientation mismatch between AE and CCX (null result)

Revision ae483fa, read-only on src/ and tests/. Probe: `logs/check/t31/probe.py` (whole es 0.25 mesh,
35,012 elements, 33,827 nodes; 88 s mesh build, probe itself seconds; no SLURM).

## 1. Orientation code paths
- T14 is not in this test's path. `tests/test_blade_ccx_parity.py:135` calls
  `PyMeshAssembler.from_model(..., list(SPAN_DIRECTION), None)` and `:272` calls
  `write_ccx_mesh(..., span_direction=(0,0,1))`. Both pass (0,0,1) explicitly, so the
  `ElementConfig.span_direction` default (`core/config.py:661`) is never consulted.
- AE: `crates/aeroelast-py/src/assembler.rs:257-296` -> `aeroelast-mesh/src/geometry.rs:37-118`
  `element_angle_offset` = angle of span projected on the element plane, measured from e1 (edge 0->1
  projected), bucketed to 0.1 deg; `assembler.rs:310-325` adds it to every ply angle.
  Fibre = projected span rotated by the ply angle.
- CCX: `core/mesh/io/writers.py:1480-1570` `_write_ccx_orientations` (RECTANGULAR, 1-axis = span (0,0,1),
  `3, angle` rotation about the shell normal), `:1590-1640` COMPOSITE section names an ORI per ply.
  CCX projects the 1-axis on the shell surface. Fibre = projected span rotated by the ply angle.
  The same definition; the bucket sets (`writers.py:88`) are only used on the non-quadratic path.

## 2. Angle statistics (AE fibre vs CCX fibre, global vectors, projected-span definition)
median 0.0004 deg, p95 0.011 deg, max 0.050 deg (= half the 0.1 deg bucket) in every span band
(z 0-5, 5-20, 20-60, 60-100, 100-118 m). No element has span parallel to its normal.
The fibre differs from the global z axis by median 1.85 deg, p95 5.3 deg, max 46.6 deg, but identically
on both sides (twist/prebend/curvature, not an AE-CCX difference).
All ply angles in all laminates are 0.0 deg, so the sign of the rotation about the normal (e3 flip) cannot matter.
Caveat: the AE vector is a NumPy port of the Rust formula read from source, not extracted from the
running extension; CCX's own internal projection is taken from its documented rule, not from an FRD.

## 3. Stiffness estimate
Not needed: no mismatch. Sign check against the hypothesis: W_AE < W_CCX means AE is the stiffer model.
A fibre rotated away from the span in AE (chordwise reading) would make AE softer in axial
(E_x falls toward E2 << E1), the opposite sign. The hypothesis is therefore also excluded by sign.

## 4. Verdict and next suspects
Null. The 9.6 -> 20.5 % growth cannot come from T14's default. Note that T14m in T21 is a
"combination" tree (T7e/T8w/T9/T14 merged), so the growth is any change in that combination
(mesh/laminate content), not specifically span_direction. Next suspects, in order:
1. Section content, not orientation: compare AE vs CCX total laminate thickness/EA per element set on
   the T14m tree (the deck writes `ply_t` from `plies_data()`, AE uses the Laminate ABD; check they
   agree set by set and that no set is missing from the deck).
2. Elements the CCX path drops or fails to integrate (`writers.py:120`; the tip sliver noted at
   `generators.py:1261`): count mesh elements vs deck elements on T14m; a CCX model missing
   load-carrying elements is softer, the observed sign for W_CCX > W_AE.
3. Quadratic (S8R) vs linear MITC4 axial response under a lumped tip load (the docstring's
   mesh table shows 3.3 % at 0.25 m; the work metric grows with local stress concentration at the tip nodes).
4. Boundary condition/load: root clamp set and tip load split across `tip` nodes on T14m
   (tip node count changes with the mesh).
Decisive next probe (needs ccx, 1 run ~ <10 min): run the T14m deck, read element count and section
thickness from the .inp and per-set EA from AE, one table.
