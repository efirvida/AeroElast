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

---
# Part 2: deck audit on T7e (f28d0ae) and T14m (ad5f3ac). Null again: the CCX deck is the same structure as our model.

Probes: `deck_probe.py`, `quality_probe.py`; logs `deck-{T7e,T14m}.log`, `quality-{T7e,T14m}.log`. Decks written exactly as
`test_blade_parity_ccx[tip_axial_tension]` (es 0.25, quadratic, span (0,0,1)) into `$SCRATCH/tmp/t31-<tree>/d.inp`. CCX not run.

| item | T7e | T14m |
|---|---|---|
| mesh nodes / elements | 32,515 / 33,673 | 33,935 / 35,120 |
| deck elements (S8R + S6) | 33,673 (33,579 + 94) | 35,120 (35,004 + 116) |
| elements dropped | 0 | 0 |
| element sets (laminate props / isotropic dict) | 9,210 (8,872 / 338) | 11,758 (11,420 / 338) |
| sets without a section in the deck; deck ELSET count != mesh | 0; 0 | 0; 0 |
| ply thickness deck vs ours (max rel) | 3.8e-7 | 3.6e-7 |
| A11 (Q11 CLT, 0 deg) deck vs ours: median / max rel | -1.5e-10 / 4.6e-7 | -1.5e-10 / 4.6e-7 |
| sum(A11*area) deck/ours | 1.000000 | 1.000000 |
| isotropic TE_FLAT dict sets: E, nu, t deck vs ours | 3.44e9, 0.3, 1 mm = same | same |
| root BC | NROOTNODES 1-6 = 0, 136 nodes (68 corner + 68 mid); ours clamps the 68 corners | same |
| tip load | NTIP DOF3 = 3.125e4 on 32 nodes (corner+mid) = 1e6; ours 62,500 on 16 corners = 1e6 | same |
| mesh z range, tip plane extent | 0..117 m, 1e-13 | same |

Ours-side EA by 20 m band (sum A11*area/20, N), T7e -> T14m: 3.552e10->3.554e10, 2.295e10->2.293e10, 2.342e10->2.342e10,
1.935e10->1.936e10, 1.017e10->1.018e10, 1.68e9->1.69e9. The axial stiffness is unchanged (< 0.6 %) between the trees.

Mesh quality (identical in both trees, so not the T7e->T14m driver): quads 33,579 / 35,004; aspect median 1.21 / 1.23,
p99 83 / 83, max 168 / 170, n(aspect>10) 656 / 670, all at z > 60 m (382 at 60-100, 274-280 at 100-117; the TE flat strips);
warp up to 0.48 of the edge length on ~61 / 63 elements; 30 / 20 elements with area < 1e-4 m2.

## What this says
1. Nothing in the deck differs from our model at section level (thickness, A11, materials, element set coverage,
   total load). The BC/load split differs only by node population (corners vs corners+midside); same total, same on both trees.
2. The case named "axial" is not axial. u_axial = F L / EA ~ 1e6 * 117 / 2.3e10 = 5 mm, but the tip response is 6.2 m
   in our model on T7e, so it is eccentric-load flexure/torsion on a prebent, twisted blade. EA is identical on both trees
   while W grows +24 % (ours) and +41 % (CCX) T7e -> T14m, so the growth is in flap/edge/torsion stiffness, and the 9.6 -> 20.5 %
   gap is the gap in those, amplified. This agrees with the S-1 GJ move (5.05 -> 10.18 %).
3. The gap sign (CCX more flexible) was not explained by orientation (part 1) or sections (above). What remains is the
   formulation: MITC4 (flat, 4-node) on warped (up to 0.48) / high-aspect (up to 170) TE-strip elements vs S8R curved
   quadratic. The set count and the TE/web changes between T7e and T14m raise the flexural share; the deck cannot tell which
   change did it. This is inference; it was not measured.

## Not found / next
No tree change in the deck explains the extra ~11 pp. The decisive discriminant needs CCX runs: run `tip_flap`, `tip_edge`,
`tip_torsion` and `tip_axial_tension` on both trees (T21's job.sbatch pattern, ~8 tasks, bigmem not needed) and compare
the components of the tip displacement (axial/flap/edge/twist) and W per case. If flap/edge/torsion gaps are each < 5 %
on T14m while the axial case is 20 %, the case is an eccentricity amplifier and should be replaced, not widened. A second
cheap AE-only test: re-solve T14m with the warped TE-flat elements (warp > 0.05) excluded from the comparison metric
or refine es to 0.125 and check whether the W gap closes (the docstring table shows axial 12.6 -> 3.3 % across es 2 -> 0.25).
