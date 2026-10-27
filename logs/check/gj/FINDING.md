# T27 FINDING: GJ gap, laminate shear (A66) vs yaml/SONATA, thick-wall Bredt

**Verdict: null result for both hypotheses.** The S-1 GJ gap is not a laminate in-plane-shear (A66) mismatch and
a thick-wall correction does not close it. At r/R 0.70 A66 contributes +0.2 pt of BD GJ, the thick-wall
open-section term +2.0 pt, and the wall-path centroid -0.2 pt. The thin Bredt is -14.8 % under BeamDyn, the
thick-wall one -12.8 %. The remaining 12.8 pt is geometry/reference (effective enclosed area), still open.

## Setup
- Tree HEAD f354972 (blade code = d5d0fd2): `BladeMesh(tests/IEA-15-240-RWT.yaml, element_size=0.25,
  skin_reference="mid")` -> `mesh_head.pkl` (`mk_mesh.py`); numad key arcs -> `arcs.pkl` (`mk_arcs.py`). No SLURM
  (mesh 90 s, each station run 5 s).
- Metric: GJ [N m^2]; reference = S-1 BeamDyn K66 (`gj_bd`), FE = S-1 `gj` from `capture-T15m.pkl` (d5d0fd2).
- "Ours": A, D of the Rust `Laminate` that `build_rust_properties` builds for the section of each cut element
  (integrated C*h, N/m). "Yaml": own windIO reader (`yamlref.py`): layers covering a nd_arc at eta, thickness and
  angle interpolated; region arcs mapped from the numad key arcs by a per-side linear fit on the cap edges; 21
  samples per region averaged. Same material law (G12 = G[0]), so a difference isolates layup/thickness/angle.
  Aeff = 1/(A^-1)66 (what Bredt uses); A66 = A[2,2]. They coincide (A16 = A26 = 0).
- Bredt: multi-cell (3 cells, 2 webs), minimum complementary energy on the actual cut of the mid-surface mesh,
  per-segment q ds / A66 (`gj_station.py`).
- Limit: the SONATA/ANBA *output* is not available per section, only BD K66 and the yaml SONATA was fed (web
  placement already verified on its PDFs, `logs/check/sonata-webs/`). "Reference A66" is therefore the ply stack
  SONATA was given. ANBA's 3D lamina law differs from the plane-stress law at percent level; not quantified.

## Station: r/R 0.70 (z = 81.9 m, mesh row 412)
S-1 FE/BD-1 (T15m): 0.44 -11.4 %, 0.50 -9.5 %, 0.56 -4.9 %, 0.62 -11.7 %, 0.685 -18.9 %, 0.746 -17.1 %,
0.81 -3.2 %. 0.70 is mid-span, inside the documented -17..-20 % dip, and where the thin Bredt itself is 15 % low.
Also run at 0.45, 0.50, 0.60, 0.80 (`gj_station_0NN.log`).

## Per-region A66, r/R 0.70 (N/m)
| section | L [m] | h ours [mm] | A66 ours | h yaml [mm] | A66 ours/yaml-1 |
|---|---|---|---|---|---|
| 01 HP_TE_REINF | 0.751 | 7.43 | 3.454e7 | 7.43 | -0.0 % |
| 02 HP_TE_PANEL | 1.098 | 17.86 | 3.555e7 | 17.84 | +0.0 % |
| 03 HP_SPAR_A/B/C | 0.094/0.693/0.099 | 72.66 | 4.402e8 | 72.62 | +0.1 % |
| 04 HP_LE_PANEL | 0.215 | 17.86 | 3.555e7 | 17.84 | +0.0 % |
| 05 HP_LE | 0.389 | 7.14 | 4.193e7 | 7.65 | +0.7 % |
| 06 LP_LE | 0.389 | 7.14 | 4.193e7 | 7.14 | +0.0 % |
| 07 LP_LE_PANEL | 0.244 | 17.86 | 3.555e7 | 16.31 | -2.5 % |
| 08 LP_SPAR_A/B/C | 0.089/0.694/0.106 | 72.66 | 4.402e8 | 72.62 | +0.1 % |
| 09 LP_TE_PANEL | 1.099 | 17.86 | 3.555e7 | 17.84 | +0.0 % |
| 10 LP_TE_REINF | 0.750 | 7.43 | 3.454e7 | 7.52 | -2.3 % |
| webs 00/01 | 0.617/0.591 | 21.88 | 2.644e7 | 21.88 | 0.0 % |

The -2.5/-2.3/+0.7 % rows are regions where the yaml thickness changes inside the region (reinforcement edge):
resolution of my region average, not a different stack. Full table with the yaml min/max: `gj_station_070.log`.

## GJ at r/R 0.70 (BD = 7.178e7 N m^2)
| quantity | GJ | vs BD |
|---|---|---|
| FE S-1 (T15m) | 5.784e7 | -19.4 % |
| Bredt thin, mesher laminate | 6.117e7 | -14.8 % |
| Bredt thin, yaml laminate | 6.133e7 | -14.6 % |
| Bredt, wall path at Q66 centroid (mean -0.2 mm) | 6.105e7 | -14.95 % |
| Bredt thick = Bredt(mid) + sum 4 D66 L (St Venant open term, +2.35 %) | 6.260e7 | -12.8 % |
| wall path on OML (upper bound) | 7.017e7 | -2.2 % |
| wall path on IML (lower bound) | 5.262e7 | -26.7 % |

Contributions, points of BD (positive closes the gap): A66 mismatch +0.2; open-section thick-wall term +2.0;
wall-path centroid -0.2 (stacks are shear-symmetric about the mid-surface); FE below its own thin Bredt -4.6;
unexplained by laminate or thick wall -12.8. Closing the Bredt gap by geometry alone needs the wall path 17.5 mm
toward the OML (1.17 x h/2, past the OML) or all-wall G x1.174. Neither is a laminate or thin-wall effect.

## Across stations (pts of BD)
| r/R | FE | Bredt thin | A66 mismatch | open term | unexplained after thick |
|---|---|---|---|---|---|
| 0.45 | -11.5 | -4.3 | +0.07 | +1.4 | -2.9 |
| 0.50 | -9.5 | -5.2 | +0.01 | +1.6 | -3.7 |
| 0.60 | -8.1 | -8.1 | +0.05 | +1.9 | -6.2 |
| 0.70 | -19.4 | -14.8 | +0.23 | +2.0 | -12.8 |
| 0.80 | -5.7 | -14.6 | -0.36 | +0.6 | -14.1 |

The Bredt-vs-BD gap grows smoothly with span (-4 -> -15 %); A66 is under 0.4 pt everywhere. FE against its own
Bredt swings from -7.2 to +9.0 pt (section/mesh dependent at es 0.25), which is the S-1 per-station scatter;
not studied here.

## Reading
1. Laminate shear: ruled out as a cause. The mesher assigns the yaml stacks within 0.1 % (caps, panels, webs).
2. Thick wall: the open term is real (+0.6..+2.4 pt; the shell carries it through D66) but small. The Q66-centroid
   path moves nothing.
3. Residual: the reference implies an effective enclosed area ~8 % above the mid-surface one (GJ ~ A^2) that
   grows outboard (inference). Next test: a direct SONATA/ANBA run of this section, which needs the tool; not
   available here.

## Files (`logs/check/gj/`)
`stations.py/.log`, `mk_mesh.py`, `mk_arcs.py` (+ `mesh_head.pkl`, `arcs.pkl`), `yamlref.py`, `gj_station.py`
(argv r/R), `gj_station_{045,050,060,070,080}.log`, `res_*.pkl`, `tab_*.pkl`, probes `api.py`, `peek.py`.
