# Auditing model parity — comparing like with like BEFORE running

Motivation: every discrepancy found so far (aligned power, yaw shape, twist,
de-loading) was measured *across* models whose configurations had never been
compared item by item.  This document is that comparison: our model vs the
reference, item by item, with a verdict per item.  No runs were launched to
produce it — it is read from the actual configuration files.

Verdict codes:

- `MATCH` — same value / same model.
- `MISMATCH` — different, and the difference is controllable.
- `CLASS` — different model class (unavoidable); the *measured* effect is
  reported so it is not confused with a defect.
- `TO CLOSE` — cannot be compared yet: we do not have the quantity on our
  side (or the definition differs and must be verified first).

---

## A. Aero side

Ours: `tests/IEA15MW/frontiersin_2025_yaw/fluid_yaw_0.yaml` +
`aeroelast/solvers/bem/engine.py` (`BEMSolver` → CCBlade, `.sources/build`
vendored copy).  Theirs: `$SCRATCH/ofruns/OpenFAST/.../IEA-15-240-RWT-Monopile_AeroDyn15.dat`
(converted to v5) + the `aerodyn_driver` deck.

| # | item | ours | AeroDyn | verdict |
|---|---|---|---|---|
| A1 | air density | 1.225 | 1.225 (`FldDens`) | MATCH |
| A2 | kinematic viscosity | 1.81e-5 | deck says `default` (≈1.4639e-5, air @15 °C); 1.81e-5 only because the driver sets it | MISMATCH vs the deck; affects Re |
| A3 | blades | 3 | 3 | MATCH |
| A4 | hub height | 150 m | 150 m | MATCH |
| A5 | hub radius | 0.0 | 1.5 m (kept from the deck) | MISMATCH (hub-loss radius) |
| A6 | precone / shaft tilt | 0 / 0 | deck −2.5 / 6; driver 0/0 (matched to ours) | MATCH in our runs, differs from the official turbine |
| A7 | wind shear | `shear_exp: 0.0` | `PLExp: 0`; deck default 0.14 | MATCH in our runs |
| A8 | wake / induction | BEM + Glauert(Buhl) | `Wake_Mod=1` (BEMT) + Glauert(Buhl) | MATCH class |
| A9 | tip / hub loss | on / on | `TipLoss=True`, `HubLoss=True` | MATCH |
| A10 | tangential induction | on (`wakerotation`) | `TanInd=True` | MATCH |
| A11 | drag in induction | CCBlade `useCd` default | `AIDrag=True`, `TIDrag=True` | TO CLOSE (confirm CCBlade default) |
| A12 | skew model | **none exists** | `Skew_Mod=1`, `SkewRedistr_Mod=1` (Glauert/Pitt/Peters), `SkewMomCorr=False` | CLASS — measured **+5.0 %** at yaw 40 |
| A13 | dynamic wake | none | `DBEMT_Mod=2`, `tau1_const=29.03 s` | CLASS — measured **≤0.3 %** |
| A14 | unsteady aero | none (quasi-steady) | `UA_Mod=3` (B-L Minnema/Pierce) | CLASS |
| A15 | tower influence/shadow/drag | none | deck `TwrPotent=1`, `TwrShadow=1`, `TwrAero=True`; driver 0/0/False | MATCH in our runs |
| A16 | nacelle drag | none | `NacelleDrag=False` | MATCH |
| A17 | blade geometry source | WindIO `tests/IEA-15-240-RWT.yaml` | `IEA-15-240-RWT_AeroDyn15_blade.dat` | MATCH provenance (same official data) |
| A18 | span discretization | 53 stations | 50 nodes | MISMATCH (interpolation) |
| A19 | sweep / curve offsets | **not used** — `BEMSolver` passes only `r, chord, theta` | `BlSwpAC`, `BlCrvAC` used | MISMATCH, TO CLOSE magnitude |
| A20 | polar / Re model | NeuralFoil at the local Re (fallback 1e7) | tabulated single design-Re table per airfoil | MISMATCH — measured ~1–3 % on aligned power |
| A21 | azimuth handling | integrated = axisymmetric average | time average (skew=0 ⇒ axisymmetric) | MATCH |
| A22 | rpm / pitch | 7.55 rpm, pitch 0, feather-positive | 7.55 rpm, pitch 0 | MATCH (unit conversion verified in code) |
| A23 | rotation / tangential sign | `tangential_direction=[-1,0,0]` | AeroDyn default | MATCH (documented in code) |

Aero summary: **no hidden convention bug.** The controllable gaps are A18,
A19, A20, A5 (and A2 in their own runs).  The model-class gaps (A12–A14)
are *measured* and small — which is why the aero is not the dominant defect.

---

## B. Structural side

Ours: WindIO `components`/`materials` (`tests/IEA-15-240-RWT.yaml`) through
`aeroelast` shell model (MITC4+/SRI, `element_size 0.25`).  Theirs: the IEA
reference — `Blade Structural Properties` sheet of `IEA-15-240-RWT_tabular.xlsx`
(the BeamDyn coordinate system: 6×6 mass `M_ij` and stiffness `K_ij` matrices
per span fraction) and the `Material Properties` sheet.

| # | item | ours | IEA reference | verdict |
|---|---|---|---|---|
| S1 | materials | WindIO `materials` (Gelcoat, steel, steel_drive, glass_uni…) | `Material Properties` sheet | MATCH (3 verified lean: Gelcoat 1235/3.44/1.323 GPa; steel; steel_drive) |
| S2 | layup definition | WindIO `components` per span | IEA layup (same provenance) | MATCH source |
| S3 | flap / edge stiffness | not computed directly (shell gives per-element A/D/B) | `K_11`, `K_22` per span (e.g. root 6.74e9 / 6.73e9 N·m²) | **TO CLOSE** |
| S4 | torsion stiffness | `K_44` equivalent only via box/ring benchmarks (S-8a, dtube) | `K_44` per span (root 1.496e11 N·m²) | **TO CLOSE (for the blade)** |
| S5 | mass per unit length | — | `M_11` per span (root 3127 kg/m) | **TO CLOSE** |
| S6 | bend–twist couplings | — | `K_14`, `K_24`, `K_46`, `K_16`, `K_26` per span | **TO CLOSE** |
| S7 | elastic axis / pitch-axis offset | our load application point (the S-8b "semichord" suspect) | BeamDyn coordinate system (EA/TC) | **TO CLOSE — prime suspect** |
| S8 | element | shell MITC4+/SRI | beam (BeamDyn) | CLASS (unavoidable) |
| S9 | mesh | 0.25 m, convergence checked | BeamDyn elements | MATCH (converged) |
| S10 | root BC | clamped | clamped | MATCH |
| S11 | **twist metric** | "section rotation" from the undeformed shell — may include cross-section distortion/warping | rigid BeamDyn section rotation | **METRIC MISMATCH — must verify before citing ×3–12** |
| S12 | load application (S-8c) | AD `Fn`/`Ft` at the local chord frame | — | MATCH |
| S13 | gravity / centrifugal | present | present | MATCH |
| S14 | damping | (not in the static comparisons) | — | TO CLOSE for dynamic cases |

Structural summary: the materials and the model provenance match, but **the
quantities that decide the twist comparison (S3–S7, and the metric S11) have
never been compared.**  The S-8c result (mid-span ×2.96, tip ×11.9 vs the
BeamDyn anchor) therefore cannot yet be cited as a structural defect: a
metric artifact (S11) or an elastic-axis offset (S7) would produce exactly
that signature.

---

## C. What to close next (all read/compute, no campaigns)

1. **S11** — read our `section_twists_deg` (in `tools/run_s7_torsion.py`) and
   state exactly what rotation it measures; compare with BeamDyn's section
   rotation definition.  If ours includes shell distortion, recompute a
   rigid-section rotation from the four section corners.
2. **S3–S6** — build our sectional stiffness/mass per span station from the
   WindIO layup and put it side by side with the BeamDyn `K_ij`/`M_ij` table
   (root, 0.25, 0.5, 0.75, 0.95).
3. **S7** — locate our elastic axis / pitch axis in the section frame and
   compare with the BeamDyn coordinate system; this is what S-8b flagged.
4. **A19** — read `BlCrvAC`/`BlSwpAC` magnitude for the IEA blade to size
   the ignored sweep/curve term.
5. **A11** — confirm CCBlade's `useCd` default.

Only after those five are closed is it legitimate to attribute the
over-twist (or the aligned offset) to a specific physical cause.

---

## D. Audit closures (read-only, 2026-09-24)

**S11 — CLOSED: the twist metric is comparable.** `section_twists_deg`
(`tools/run_s7_torsion.py:112`) is a *least-squares rigid rotation*
`θ = Σ(r×u)/Σr²` over each section's nodes, not a raw chord measure.  It is
the same construction as S-1 and equivalent in kind to BeamDyn's rigid
section rotation: shell distortion does not enter the fitted rotation.
⇒ The S-8c ratio (mid-span ×2.96, tip ×11.9) is a real twist difference,
not a metric artifact.

**S4 — CLOSED, and it rules torsion OUT as the cause.** `test_iea15mw_s7_torsion.py`
documents: our shell GJ is **+17–19 % stiffer** than the BeamDyn tables
(twist ratio ≈ 0.84 under a pure moment).  So under equal torsion we twist
*less*, not more.  The over-twist under aerodynamic loads cannot come from
the torsional stiffness.

**S/load-point — the prime suspect, now located.** `apply_ad_loads`
(`tools/run_s8c_twist_sign.py:35`) applies the AD `Fn`/`Ft` as nodal forces
on **every shell node** (`u_x, u_y`), so the resultant of each station acts
at the **node centroid of the cross-section**.  AeroDyn/BeamDyn carry the
resultant at the **aerodynamic center** plus the airfoil **pitching moment
`Cm`**, which the S-8c does not apply at all.  Two distinct consequences:

1. the chordwise application point differs (node centroid vs AC / pitch
   axis) — this is exactly what S-8b flagged;
2. the aerodynamic pitching moment term is missing entirely.

Note the sign test: omitting `Cm` (nose-down, negative for these airfoils)
removes moment, so it cannot explain an *excess* twist.  That leaves the
application point and/or our section's **bend–twist coupling**.

**Consequence for the causal chain:** the earlier claim "the structure
over-twists, therefore the coupling over-de-loads" is still open.  What is
closed: it is not GJ, and it is not the twist metric.

## E. Remaining items, in priority order (all read/compute)

1. **S6 — our section bend–twist coupling** (equivalent of `K_14`, `K_16`,
   `K_46`) built from the WindIO layup, side by side with the BeamDyn table.
   This is now the top suspect.
2. **S7 — the torsional/elastic axis position** in the section frame vs the
   BeamDyn coordinate system, and the resultant's arm from it.
3. **S6b — whether the `Cm` term exists at all** in `extract_ad_loads` /
   `apply_ad_loads` (its absence must at least be declared in the docs).
4. **A19 — `BlCrvAC`/`BlSwpAC` magnitude** for the IEA blade (sweep/curve we
   ignore).
5. **A11 — CCBlade `useCd` default** confirmation.

---

## F. ERRATUM — A20 was wrong, and Peldaño 2 must be re-read

**Correction, with evidence.** `load_blade_aero` (`aerodynamics.py:583-588`)
states: *"parses the ``outer_shape_bem`` section for chord, twist, pitch axis
and reference axis, and the ``airfoils`` section for polar data. When an
airfoil entry has no ``polars`` key, NeuralFoil is used as a fallback"*.
The WindIO file **does** carry `polars:` blocks (`tests/IEA-15-240-RWT.yaml`
lines 565-715): circular `re 3.0e6`, SNL-FFA-W3-500 `re 8.1e6`, and the
FFA-W3-* family at `re 1.0e7`.  Therefore:

- **ours = the IEA tabulated polars**, not NeuralFoil.  `neuralfoil_model`
  in the fluid yaml only affects the fallback path.
- My table entry A20 ("ours = NeuralFoil at the local Re") is **wrong**.
- Consequently the **Peldaño 2 conclusion is withdrawn**: what that script
  actually compared was the *WindIO-embedded IEA tables* against the
  *AeroDyn polar files*, and it evaluated ours at `re=1e7` while every
  AeroDyn file is tabulated at **`re=3.0e6`** (verified for files 00, 04,
  16, 28, 41, 49 — all 3.0e6).  The +0.3…+4.3 % CL deltas are a
  **Reynolds-selection artifact of that comparison**, not a polar model
  difference.

**Corrected A20**: the polar *source* is the same IEA data on both sides;
the live mismatch is the **Re selection** — ours uses the per-airfoil table
(3e6 / 8.1e6 / 1e7) with CCBlade's local-Re lookup, AeroDyn uses the single
**3e6** table for every station.  Magnitude on aligned power: small — the
AeroDyn run (15.815 MW) lands within ~0.4 % of the IEA `Rotor Performance`
sheet, so the Re choice does not carry the +4.4 % gap.

**What this does not change:** the structural closures (S11, S4) and the
load-application finding are on a different axis and stand.

**What is now unattributed:** the R1-vs-R2 aligned gap (+4.4 %) is still
open; its remaining candidates are A18 (53 vs 50 stations), A19 (ignored
sweep/curve), A5 (hub radius 0 vs 1.5) and this Re selection — all
checkable without new campaigns.

---

## G. A18/A5 CLOSED — the blades are geometrically different

Three representations of the same turbine, compared at equal `r/R`
(non-dimensional, so the datum cancels):

| r/R | chord OFF | chord AD | chord OURS | OURS−OFF | AD−OFF | twist OFF | AD | OURS | OURS−OFF |
|---|---|---|---|---|---|---|---|---|---|
| 0.15 | 5.648 | 5.607 | 5.535 | −1.99 % | −0.73 % | 11.03 | 11.58 | 12.47 | **+1.44°** |
| 0.30 | 5.368 | 5.433 | 5.531 | **+3.04 %** | +1.21 % | 5.52 | 5.78 | 6.21 | **+0.69°** |
| 0.50 | 4.155 | 4.188 | 4.241 | +2.07 % | +0.81 % | 1.69 | 1.78 | 1.91 | +0.22° |
| 0.70 | 3.221 | 3.238 | 3.263 | +1.31 % | +0.54 % | −0.49 | −0.44 | −0.37 | +0.12° |
| 0.85 | 2.525 | 2.535 | 2.545 | +0.79 % | +0.37 % | −2.17 | −2.17 | −2.17 | 0° |
| 1.00 | 0.500 | 0.500 | 0.500 | 0 % | 0 % | −1.24 | −1.24 | −1.24 | 0° |

Rotor radius: **official blade coordinate 117.0 m**, **AeroDyn 118.5 m** (BlSpn
+ 1.5 hub), **ours 121.12 m** — our rotor is **2.2 % larger** than the deck's.

Sources: official = the `Blade Geometry` sheet of
`IEA-15-240-RWT_tabular.xlsx` (Blade Span [r/R], Rotor Coordinate, Chord,
Twist, Prebend, Sweep); AD = `IEA-15-240-RWT_AeroDyn15_blade.dat`;
ours = `tests/IEA-15-240-RWT.yaml` through `load_blade_aero`.

**Reading.** The AeroDyn blade tracks the official sheet to 0.2–1.2 %.
**Our blade does not**: chord +0.1…+3.0 %, twist up to **+1.44°** inboard,
radius +2.2 % vs the deck / +3.5 % vs the official blade coordinate.  These
are first-order power drivers (chord → torque; twist → local AoA; radius →
TSA and tip loss), so the R1-vs-R2 aligned gap (+4.4 %) is **substantially a
geometry mismatch**, not an aero-model difference.

Our YAML is an official IEA WindIO file ("with taped chord tip design" is in
its own `name`), so both sides are defensible *representations* of the IEA
15 MW — but they are not the same blade, and comparing them without saying
so was itself a parity error.

**Consequences to decide (not to guess):**

1. For the aero ladder to mean anything, both sides must use the **same
   blade definition**.  Cheapest path: feed our (r, chord, twist) into the
   AeroDyn blade file, or the AeroDyn blade into our loader.
2. The campaigns currently run on our YAML geometry (R = 121.12 m, chord
   +3 %).  Whatever the article reports, it must state which blade
   representation is the reference and why.
3. `BlCrvAC` max = **3.999 m** and `BlSwpAC` max = 0.435 m (A19): the blade
   has a prebend we ignore entirely.  For a 117 m blade that is ~3.4 % of
   the length — not negligible for the local inflow and for the structural
   pre-bend coupling.

---

## H. S7 and S6b CLOSED — three parity errors, the ×3–12 chain is withdrawn

**S7 — the pitch axis: present in the data, carried to the stations, not used
by the aero.**
`tests/IEA-15-240-RWT.yaml:24` defines `pitch_axis` as a grid/values pair;
`aerodynamics.py:672-673` reads it and `:769` hands it to each station;
`fsi_participant.py:866` re-supplies it when rebuilding stations on the
deformed geometry.  But `BEMSolver` builds CCBlade from `(r, chord, theta)`
only — **CCBlade has no pitch-axis argument**.  So on the aero side the pitch
axis is carried and unused; the equivalent physical statement is made by the
force projection onto the mesh.  `aerodynamics.py:555` also shows one branch
setting `pitch_axis` from the *aerodynamic center* — a conflation worth
keeping in mind, though not on the path used by the campaign.

**S6b — the pitching moment is computed but never applied offline.**
`engine.py:339` computes `Mp = q_c · chord · cm` per station, so the coupled
solver does carry the aerodynamic pitching moment into the projection.  But
the S-5/S-8c offline chain drops it: the reference CSV
(`tests/IEA15MW/reference/s5_ad_blade1_loads.csv`) has only
`r_m, fn_Nm, ft_Nm, twist_deg`, and `extract_ad_loads` (`run_s5_oneway.py:113`)
documents exactly those four returns — **the AeroDyn node pitching moment is
never read.**

**Why this matters, with the sign.**  The FFA-W3 sections carry a nose-down
(negative) `Cm` at rated.  Omitting it removes a twist-*reducing* moment, so
the applied load case is *more nose-up* than reality — which is the direction
of the excess I measured.  Combined with loads applied at the node centroid
rather than the aerodynamic center (section D), the S-8c is **not a faithful
reproduction of the FSI load case**.

### Withdrawn

> "The structure over-twists ×2.96 (mid-span) / ×11.9 (tip), therefore the
> coupling over-de-loads."

**Retracted.**  The S-8c omits the pitching moment and misplaces the
resultant; both are first-order for a torsion measurement.  The twist
comparison must be rebuilt with the full load case (Fn, Ft **and** the node
Mm, applied at the aerodynamic center) before any structural conclusion is
drawn.

### What still stands, from measurements that have no such flaw

- **FSI over-de-loading ×2.0–2.7** (País 3): measured on actual coupled runs
  against the literature — no missing term involved.
- **GJ +17–19 % stiffer than BeamDyn** (S-7 test): a direct torsion
  benchmark.  Note the direction — a *stiffer* blade twists *less*, so a soft
  structure cannot be the cause of the over-de-loading either.
- Every aero parity item listed in sections A/F/G.

### Meta-lesson (three for three)

Three "mismatches" in this audit were comparisons without parity: the polar
model (F), the blade geometry (G), and the load case (H).  Each one would
have sent us to fix something that was not broken.  The parity rule is now:
no comparison is reported without first stating, for both sides, the blade
geometry, the polar table and Re, and the complete load case.

### Next diagnostic that is *not* another offline shortcut

Instrument the FSI participant itself (per converged window: delta-twist from
`_compute_deformed_geometry`, local Re, CL, alpha, and the applied `Mp`) and
compare each against the rigid-window values.  That localises the
de-loading *inside* the coupling with the real load path, instead of
inferring it from a hand-built static case.
