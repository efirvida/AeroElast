# T24 open mass findings — located causes (no code changed)

Revision: HEAD d8a1c3a, code d5d0fd2 (mesh-affecting code unchanged since). Official `tests/IEA-15-240-RWT.yaml`, es 0.25.
Scratch (this dir): `bins-head.json` / `bins-oml.json` (`mass_bins_ref.py <yaml> <out> mid|oml`, T15 method), `yaml_mass_layers.py`,
`stations.py`, `materials_by_station.py`, `span_integral.py/.csv/.log`, `layup_vs_laminate.py/.json`, `audit/` (`tools/audit_blade_input.py`
at HEAD, OML), `ref_compare.txt`. Mesh pickle: `$SCRATCH/tmp/t24-mesh.pkl` (outside the repo).

## Baseline reproduced (skin_reference="mid", the anchors' setting)
M 67,150.0 / skin 64,648.8 / webs 2,500.6 kg; yaml skin 64,842.2, webs 2,627.6 (T15 REPORT identical). TE_FLAT glass_uni 157.0 + 209.7 = 366.7 kg;
glass_uni all 5,889.0 vs yaml 5,644.5 (+4.3 %). Triax -215.4, carbon -213.9, webs -127.0 (-4.83 %).

## Key discovery: the T15 comparison mixes two reference surfaces
The yaml integral uses **outer-surface** arcs and straight outer-point web heights. The mesh in `mass_bins.py`/V-03/S-x is
`skin_reference="mid"` (`core/mesh/generators.py:1181` -> `core/mesh/skin_reference.py:8`, insets every outer node by the mean half-thickness).
The same script with `oml` (`ref_compare.txt`):

| kg | yaml | OML mesh | OML vs yaml | mid mesh | mid vs yaml | mid - OML |
|---|---|---|---|---|---|---|
| triax | 40,475.1 | 40,801.4 | +0.81 % | 40,259.7 | -0.53 % | -541.7 |
| carbon | 14,895.4 | 14,876.0 | -0.13 % | 14,681.5 | -1.44 % | -194.5 |
| webs | 2,627.6 | 2,635.5 | +0.30 % | 2,500.6 | -4.83 % | -134.9 |
| glass_uni TE_FLAT | 0 (not in yaml arcs) | 366.7 | | 366.7 | | 0 (TE_FLAT stays on OML, `skin_reference.py:61`) |

Independent check (planar cuts of the OML pickle, 1 m midpoint rule, `span_integral.log`): triax +0.59 %, carbon -0.22 %, webs +0.19 %;
mesh web height equals the yaml outer-point height to 0.3 % at z=20/40/60 (`stations.py`); web thickness, attach nd_arc and ply order pass F4;
set layup == solver laminate (`layup_vs_laminate.py`: ratio 1.0000 on 10,250 skin + 1,170 web sets), so no density/unit mismatch.

## (a) TE_FLAT glass_uni 367 kg — translation defect (mesher region mapping)
- yaml: `TE_reinforcement_SS` start_nd_arc fixed `TE` (0), width 0.75 m, t 30 mm (z>~20 m); `_PS` end_nd_arc fixed `TE` (1). nd_arc excludes the flat-back face
  (audit convention). Yaml arcs at z=20/40/60: [0,0.0561]/[0,0.0676]/[0,0.0856] = 0.750 m.
- mesher: `yaml_to_blade.py:449-450,454,460` assigns the band extents `["d","te"]`; `Component.find_region_extents` (`component.py:98`) maps it to key
  segments 0..2 = `TE_FLAT`(te-e) + `TE_REINF`(e-d); `stackdb.py:91,107` adds the ply to every segment in range. The band is not "carried by an arc mapping"; the extent
  label `te` simply includes the flat face. Mesh at z=40: HP_TE_REINF cut 0.749 m, LP 0.750 m (= yaml 0.750), plus HP/LP_TE_FLAT cut 0.060 m each with the full 30.00 mm glass_uni
  (z=20: 0.056/0.126 m, 30.43 mm; z=60: 0.029 m each, 30.00 mm). The flat face is real mesh geometry (x const, y varies).
- Impact: 366.7 kg = 0.55 % of blade mass; glass_uni +4.3 % -> -2.2 % (-121 kg) without it (remaining shortfall is TE_REINF, -1.0/-3.4 %, short at z<30: reinf cut 0.683/0.616 m at z=20 vs 0.75).
  Stiffness (ring integral E1 t ds, `materials_by_station.py`): EA +0.5-1.1 %, **EI_edge +3.1 % (z 15), 5.7 % (20), 4.7 % (30), 6.9 % (40), 3.6 % (60), 0.6 % (80)** — the flat face sits 3 m from the centroid. EI_edge is currently
  -1.9 % vs BeamDyn, so removing it would push edge stiffness lower (to about -6 %): any fix must be re-judged on S-4 (1st edge already -5.6 %, rejected at 5 %), V-02, S-2 together, not on mass alone.
- Classification: **translation defect (mesher-side region choice inherited from NuMAD's `te` extent); yaml does not place it.** Whether the physical blade carries glass on the flat back: unresolved (reference decks not inspected).
- Candidate fix (not implemented): for the TE reinforcement/shell-band extents use the `e` label (`["d","e"]`) so TE_FLAT gets no band ply, or drop `glass_uni` from `*_TE_FLAT` stacks; the 5 % edge anchor decides. Mesh change -> restamp all blade rows.

## (b) triax -215 kg, carbon -214 kg — reference-surface (mid-surface arc shortening) + an open OML root surplus
- Not arc truncation, not density/unit, not G10 span interpolation: on the OML mesh carbon is -0.13 % and per-station line mass matches the yaml within ~0.5 % for z>=20 (`materials_by_station.py`).
- Where: carbon deficit is spread z 10-100 (mid, per 10 m: -18, -29, -28, -47, -10, -53, +1, -28); OML carbon is within ±25 kg per 10 m. Triax: mid -98/-97 kg at z 0-20, then ~0 outboard.
  Cause: the inset moves each outer node by the half thickness of its skin laminate, so the mid-surface arc of a layer is shorter than the yaml outer arc (carbon-bearing spar caps are thick and curved; root is a thick cylinder: perimeter 16.11 m mid vs 16.36 yaml open polyline at z=2).
- Residual: the OML mesh is **+265 kg (0-10 m) / +38 kg (10-20 m)** over the yaml for triax (cut-based: +1.1 to +2.4 % at z 3-10, +45/+122/+79 kg per 5 m bin at 0-15). That surplus masks about 60 % of the mid-surface loss (net mid -215). Cause not located: it tracks the steep root taper (candidate: per-bay thickness sampling of a steep triax taper, or the root arc mapping); needs a bay-vs-linear thickness check at z 3-15.
- Impact: mid shortening -542 kg triax / -195 kg carbon; stiffness effect on carbon is in EA/EI by the same -1.3 % (a real shell-reference effect, shared by the ring integrals that F5 passes at 0.30/0.60).
- Classification: **yaml-faithful (reference-surface effect; the yaml outer-arc integral overstates the layer mass)** for the mid-vs-yaml part; **root OML surplus unresolved** (+0.8 % triax, 303 kg).
- Candidate: none needed for the mid part; compare against a mid-surface yaml reference (arc at ply mid-depth) so the audit stops reporting it as a deficit; investigate the root bay sampling separately.

## (c) webs -4.8 % — yaml-faithful (reference surface), not a web defect
- OML: +0.30 % (2,635.5 vs 2,627.6); mid: 2,500.6 (-4.83 %); biax -4.87 %, foam -4.78 % (same ratio -> geometry, not layer thickness). Web heights mesh/yaml at z=20: 2.5957+2.6891 vs 2.6041+2.6806 m; z=40: 1.578+1.637 vs 1.5825+1.6407 m; z=60: matches (OML). Thickness (F4) and span limits match.
- Cause: the web is sized to the yaml outer-point height; the inset moves its end nodes by the half skin thickness at the spar caps -> web height shorter by about (t_HP+t_LP)/2 (-134.9 kg, concentrated z 10-70 where caps are thick). Web span limits are not involved (the mesh has webs from z≈10, where the yaml web thickness is zero).
- Classification: **yaml-faithful (mid-surface by design, `beamdyn-arbiter-midsurface.md` T4)**. Not a translation defect; the yaml web height is an outer-point distance, an upper bound for the physical web (bonded to the inner skin faces).
- Impact: -127 kg (0.19 % of blade); webs carry shear/torsion (GJ gap): the 4.8 % is a shorter web in the mid model, the same stiffness choice the anchors were calibrated on.
- Candidate: none (mid is the anchor reference); optional: report web mass vs a mid-surface yaml height (h - (t_HP+t_LP)/2) in the audit.

## Verdict
| item | reproduced | cause | impact | class | candidate fix |
|---|---|---|---|---|---|
| a | 366.7 kg glass_uni in TE_FLAT | `yaml_to_blade.py:449-460` extent `te` -> `component.py:98`, `stackdb.py:91` | +0.55 % mass, EI_edge +3-7 % z15-60 | translation defect | extent `e` / drop glass from TE_FLAT; re-judge on S-4/V-02/S-2 |
| b | triax -215, carbon -214 kg | mid inset `skin_reference.py:8` (-542/-195); OML root surplus +303 triax | -0.5 / -1.4 % | yaml-faithful + unresolved root | mid-depth yaml reference; root bay check |
| c | webs -127 kg (-4.83 %) | mid inset shortens web by (t_HP+t_LP)/2 | -0.19 % blade | yaml-faithful | none |

Caveats: element-centroid 1 m bins (mass_bins) vs 1 m planar cuts (audit slice) agree on OML totals within 0.5 %; es 0.25 single size; EI numbers are ring-only E1 t ds (webs excluded), orientation ignored (fibre along span).
The audit tool (`tools/audit_blade_input.py`, `Blade(...)` default `oml`) and the anchors (`mid`) use different reference surfaces; the audit's F4/F5 passes say nothing about the mid mesh.
