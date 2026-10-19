# rotating-local-tip-modes

Why do low-beam-fraction (local) modes appear at 7.56 rpm but not parked, with the residual at
z/L 0.94-0.97? IEA-15 official blade, es 0.25, oml and mid. Follows `odd/tasks/s4-mode-classifier.md`.

## Specs

- S1. "Abrir los `.npz`, ubicar los nodos con residuo y ver si coinciden con el cambio de tamaño de malla en la punta o con un borde libre."
- S2. "Recuperar las tensiones de membrana centrífugas en esa zona y compararlas con la referencia analítica (σ ≈ ρω²∫r dA / A). Si da compresión o un pico donde debería haber tracción suave, el problema está en K_G."
- S3. "Repetir el modal con la rotación prendida y K_G apagado, y con K_SP apagado, para aislar cuál de los dos genera los modos."
- S4. "Mirar la tendencia al refinar la malla: 0.25 contra 0.15 m. Un artefacto numérico se mueve con la malla; un modo físico converge."
- S6. "si armalas y ponlas a correr, recuerda que tambien tenemos openseespy que podria ser otro arbitro tambien si fuera necesario" (context: CalculiX reference proposed in L7: static CENTRIF, then *BUCKLE and *FREQUENCY PERTURBATION on the same shell mesh).
- S7. "vamos con la vlalidacion analitica en otra sesion" (context: L15 proposal, T7 i-iv).
- S5. "El rojo de S-4 en el 1er modo de borde (5.9–6.1 % contra un límite de 5 %) es un hallazgo de precisión. No se toca la banda."

## Tasks

- T1 [S1] locate local-mode residual (elements, section, skin/web, tip refinement stations) - inline - DONE (L7) - commit: (tools, this handoff)
- T2 [S2] centrifugal prestress at 7.56 rpm: per-element membrane resultants N vs the section tension T(r) - inline + sbatch - DONE (L7) - commit: -
- T3 [S3] rated modal with K_G off, and with K_SP off - inline + sbatch - DONE (L4) - commit: -
- T4 [S4] es 0.25 vs 0.15 local-mode frequencies and beam fractions - sbatch - DONE (L6) - commit: -
- T6 [S6] independent shell arbiter for the tip local modes: our linear buckling factor + CalculiX S8R CENTRIF -> *BUCKLE and *FREQUENCY PERTURBATION at 7.56 rpm, same mesh es 0.25 - worker + sbatch - commit: -
- T7 [S7] analytic plate-buckling anchor for MITC4 K_G (L15 plan i-iv) - NEW SESSION - commit: -
- T5 [S5] record finding (closures doc, S-4 docstring if a number moves), close - inline - commit: -

## Log

- L1 (2026-10-18, maintainer): "vamos con ese plan" - the plan is the four experiments S1-S4 plus S5,
  proposed after `mem_search "s4-mode-classifier"`. Commits on `main`, not pushed (maintainer practice).
- L2 (evidence, before T1): K_SP cannot be the cause. Parked, no local mode below 3.8 Hz
  (omega^2 >= 580 rad^2/s^2); rated oml local at 1.5325 Hz (omega^2 = 92.7). K_SP shifts omega^2 by
  at most -Omega^2 = -0.627 rad^2/s^2. A drop of that size needs a compressive K_G locally (real or
  spurious). K_G path: `assemble_geometric_stiffness` -> static `K u = f_cf` -> `compute_element_stress`
  (membrane, centroid, `cm_raw = A/h`, so sigma = N/h, B*kappa coupling dropped) -> `assemble_geometric_k`
  (x h). The mesh adds 2 interpolated tip stations (`_refine_high_gradient_sections`).
- L3 (tools, uncommitted): `tools/diagnose_rotating_prestress.py` (prestress N, section force vs
  exact T(z), per-element Rayleigh split of the stored rated modes, K_G COO split per element);
  `assemble_rotating(..., include_kg, include_ksp)` + `diagnose_mode_beam_projection.py --without`.
  SLURM: T1/T2 11618436 (oml) 11618437 (mid) -> logs/check/tipmodes/prestress-*; T3 11618438
  (oml no kg) 11618439 (oml no ksp) 11618440 (mid no kg); T4 11618441 (oml es 0.15, rpm 0 7.56).
- L4 (2026-10-19) T3 RESULT, es 0.25, 7.56 rpm (`logs/check/tipmodes/analyse-T3*.txt`):
  * K_G is the sole cause. oml without K_SP reproduces the full rotating spectrum to <=0.1 %
    (local 1.5325/2.3636/2.8491/3.4524 Hz, beam frac 0.15-0.21, residual peak z/L 0.94-0.97).
    oml without K_G reproduces the parked spectrum (1st flap 0.5524 vs parked 0.5534; first
    local 3.834 Hz = the parked local, residual peak z/L 0.84, not the tip) and every tip local mode
    below 3.8 Hz disappears. mid behaves the same (rated locals 2.168/2.918/2.963/3.610 Hz -> none
    below 3.80 Hz without K_G).
  * The local frequencies move strongly with skin_reference (oml first local 1.53 Hz, mid 2.17 Hz),
    while the beam modes move <=2 %: consistent with the open B*kappa hypothesis, not proof of it.
  * Caveat: removing K_G also removes the legitimate stiffening (1st flap 0.574 -> 0.552 Hz), so
    T3 identifies the term, not the defect inside it. That is T1/T2.
- L5 (2026-10-19) T1/T2 jobs 11618436/37 died after the prestress solve (residual 5e-10):
  `asm.elements` now yields `_ElemProxy` (.id, .node_ids) - no `.nodes`; `getattr(el, "thickness")`
  would also have silently given h = NaN. Fixed in the tool (MeshElement by id, h = layup sum,
  SystemExit on missing h); resubmitted as 11618652 (oml) / 11618653 (mid).
- L6 (2026-10-19) T4 RESULT, oml 7.56 rpm, es 0.25 -> 0.15 (`logs/check/tipmodes/analyse-T4.txt`):
  * Beam modes converge: 1st flap 0.5742 -> 0.5746 (+0.06 %), 1st edge 0.7044 -> 0.7044,
    2nd flap 1.6521 -> 1.6430 (-0.55 %), 2nd edge 2.0947 -> 2.1051 (+0.50 %).
  * Rated tip locals do NOT converge: 1.533/2.364/2.849/3.452 Hz (4 below 3.5 Hz) become
    2.053/3.323 Hz (2 below 3.5 Hz); the lowest one moves +34 %. Parked locals move only
    3.834 -> 3.692 Hz (-3.7 %). With T3 this reads as a mesh-dependent K_G artefact, not physics,
    but two meshes are not a trend: an es 0.125 point would close it.
  * Hazard for S-4: at es 0.15 the 2.053 Hz local sits 2.5 % below 2nd edge, whose beam fraction
    drops 0.993 -> 0.659 by mixing; any classifier keyed on beam fraction can mislabel 2nd edge.

- L7 (2026-10-19) T1/T2 RESULT, es 0.25, 7.56 rpm. Jobs 11618652/53 were moved to sequana_cpu_bigmem
  and COMPLETED (prestress residual 5e-10). Analyse logs: `logs/check/tipmodes/analyse-{oml,mid}-025.log`.
  Sanity: (a) `rq == w2` in every mode, both refs: PASS. (b) The per-bin `T_fe/T_ref` scatters 0.74-1.25
  (std 0.137, mean 0.9985), with the same pattern in oml and mid. This is a tool artifact: a 0.585 m bin
  catches 2 or 3 rows of 0.25 m elements by centroid. Window integrals int T_fe / int T_ref: 0.9999 /
  1.0035 / 1.0012 / 0.9997 / 0.9963 / 0.9927 over z/L 0-0.9, and 0.982 / 0.987 / 1.005 at the tip bands:
  PASS. The recovered N carries the exact axial resultant, so neither 1/h nor a lost tension is involved.
  Findings: (1) The local modes (beam frac 0.12-0.44) are near-cancellation modes: K_G removes 90-95 %
  of the elastic Rayleigh stiffness (oml 1.53 Hz: K=1756, KG=-1664; mid 2.17 Hz: K=1871, KG=-1686).
  K_SP < 0.1. The Rayleigh buckling factor of these shapes is K/|KG| ~ 1.05. (2) The softening is
  distributed, not slivers: the top 12 elements carry 13-27 % of the negative K_G energy. They sit in
  `09_45..09_47_LP_TE_PANEL` (h 2.7 mm, area 0.017-0.024 m2) at z/L 0.92-0.97, plus `10_38_LP_TE_REINF`
  at z/L 0.79. (3) Compression is real in the recovered field. Nmin < 0 in 79-86 % of the tip elements.
  Nmin p01 is -4.3 kN/m at z/L [0.94,0.97), and Nzz p01 ~ Nmin p01, so the compression is spanwise.
  Inboard Nzz p01 is -2.1e5 N/m: centrifugal-induced bending on top of the tension. (4) oml and mid
  have nearly the same prestress (band stats within 1-4 %) and the same K, KG magnitudes (within 6 %).
  The 2x difference in net stiffness (92.7 vs 185.5) is the near-cancellation amplifying the K
  difference. This also explains the T4 mesh sensitivity. It does not test B*kappa locally: the omitted
  B*kappa term integrates to ~0 over a section, and its local size is unmeasured.
  Open question: is the TE panel near the tip really at ~95 % of local buckling under centrifugal load
  alone? ED and BD cannot answer it (no cross-section deformation; see Engram
  `odd/rotating-local-tip-modes/reference-kg-comparison`). Candidate independent reference: CalculiX on
  the same mesh, static CENTRIF, then *BUCKLE (buckling factor vs our K/|KG| ~ 1.05) and *FREQUENCY
  PERTURBATION.
- L8 (2026-10-19, maintainer) S6 verbatim -> T6. OpenSeesPy NOT attempted: documented as not viable for the
  blade (no +-45 deg plies in this build; GLIBC 2.28 needs singularity; `odd/tasks/validation-l3-b3-audit.md`
  L82/304, `odd/tasks/beamdyn-arbiter-midsurface.md` T10); reported to the maintainer.
  PREDECLARED before any T6 run (decision criteria, not pass bands):
  (a) ours: lowest positive linear buckling factor lam of K phi = lam (-K_G) phi from the 7.56 rpm
      prestress, oml and mid. Consistency with L7: expect lam ~ 1.0-1.1 (Rayleigh K/|KG| ~ 1.05 is an
      upper bound for the lowest mode).
  (b) CalculiX S8R, same es 0.25 mesh, the T10 realisations: oml laminate centred (implementation check
      vs aeroelast oml) and oml *SHELL SECTION OFFSET (vs aeroelast mid). Prestress check first: section
      tension within 2 % of T_ref (window integrals as L7) else the CENTRIF setup is wrong.
  (c) PHYSICAL if CalculiX gives lam_ccx <= 1.5 with the buckling shape in the tip TE panels (z/L >= 0.9)
      AND rotating local modes below 3.8 Hz in *FREQUENCY PERTURBATION. ARTIFACT if lam_ccx >= 3 and no
      local mode below 3.8 Hz. Anything between is reported as unresolved, with the numbers. Near
      cancellation (L7-4) means local frequencies are not compared to a tolerance; lam is the metric.
- L9 (2026-10-19) T6 built (worker mvf92bb0-1-xojt + parent fix) and SUBMITTED on sequana_cpu_bigmem:
  11618694 buckle-oml, 11618695 buckle-mid (`diagnose_rotating_prestress.py buckle`, (-K_G)phi = mu K phi),
  11618696 t6-ccx-oml, 11618697 t6-ccx-offset (`tools/t6_ccx_centrifugal_buckling.py`, out
  `logs/check/tipmodes/t6-ccx-{oml,offset}/result.json`). CCX conventions from the tiny smoke
  (`logs/check/tipmodes/t6-syntax-smoke.py`, syntax evidence only): CENTRIF axis +y, centre 0, no hub offset
  (= diagnose tool :89); BUCKLE factors are total (do not add 1); BUCKLE .frd = static + time-0 reference +
  modes (take the last n blocks); CCX also reports NEGATIVE factors (load reversal, ignore; deck asks 20 so
  they do not crowd out the positives); OUTPUT=2D must be set in the static step. The CCX tool withholds the
  eigen runs if the L8(b) prestress windows miss 2 %.
- L10 (2026-10-19) Buckles moved to sequana_cpu_dev (maintainer: "esas corridas que son cortas y de bajos
  recursos ponlas en sequana_cpu_dev siempre que necesiten menos de 20min"; dev allows 1 submitted job).
  Dev job 11618699 oml took 7m10s but is INVALID: the regular-mode eigsh ((-K_G)phi = mu K phi, K metric)
  returned unconverged pairs (relative residual 0.999, lam 0.006-0.010, all at one node z/L 0.80);
  file kept as `buckle-oml-025-INVALID-unconverged.npz`. Fixed: ARPACK buckling mode, shift-invert
  about sigma (default 0.5) with (K - sigma Gf) factored by MUMPS; verified on a synthetic SPD K /
  indefinite G pair (factors equal to the dense reference, residual 3e-15). Re-run as one dev job per
  reference, in sequence (oml, then mid).
- L11 (2026-10-19) Maintainer: "quita el calculix de bigmeme y encolalo en dev" -> CCX bigmem jobs
  11618696/97 cancelled; `t6_ccx_centrifugal_buckling.py --stage {static,buckle,frequency}` (one stage per
  20-min dev job; eigen stages refuse without a valid static stage). The symmetric buckling-mode eigsh
  still failed (11618724 oml 2m27s, 11618727 mid 1m23s: residual 0.92-0.99, lam clustered at sigma),
  EXCEPT one oml pair: lam = 0.5177, residual 5.7e-10, peak z/L 0.864 in 09_42_LP_TE_PANEL. A pair with
  that residual is a genuine eigenpair: if confirmed, K + K_G is INDEFINITE at 7.56 rpm (panel above its
  linear critical prestress). Not claimed until the fixed solver confirms it. Fix: general Arnoldi on
  (K - sigma Gf)^-1 Gf, lam = sigma + 1/nu, plus a printed asymmetry of K and -K_G (hypothesis: Lanczos
  failed on non-symmetric matrices); verified on a slightly non-symmetric synthetic pair (residual 1e-15).
  One driver (job-6) runs the dev queue in sequence: buckle oml, mid, then CCX oml/offset x 3 stages.
- L12 (2026-10-19) T6(a) RESULT, our linear buckling, es 0.25, 7.56 rpm prestress (Arnoldi; 11618728 oml
  2m35s, 11618729 mid 1m06s; logs `buckle-{oml,mid}-025-dev-1161872{8,9}.log`). Asymmetry of K and -K_G is
  1e-17, which REFUTES the L11 hypothesis; the cause of the Lanczos failure is unexplained, and Arnoldi is
  used regardless. All 10 pairs per reference converged (residual 1e-11..2e-10).
  oml lam = 0.518 0.576 0.644 0.734 0.843 0.874 0.955 1.049 1.114 1.184;
  mid lam = 0.541 0.599 0.662 0.761 0.876 0.883 0.922 0.987 1.101 1.170.
  The modes are LP_TE_PANEL 09_42..09_47 at z/L 0.86-0.97, beam fraction 0.08-0.14, except mid #6
  (HP/LP_TE_FLAT, z/L 0.80). L8(a) expected lam ~ 1.0-1.1: the measured lam_min ~ 0.52-0.54 lies BELOW it
  (still consistent with the Rayleigh upper bound). Consequence: by our model, the rated centrifugal
  prestress is ~1.9x the linear critical load of the tip TE panels, so K + K_G is INDEFINITE (7 oml / 8 mid
  factors < 1). A rotating modal that searches only above a positive shift misses these negative-
  stiffness modes; this is an open S-4 hazard to check. CCX stage static failed in the parser (S6 web
  triangles: 6 points per ply); fixed and the 6 CCX stages re-queued on dev.
- L13 (2026-10-19) T6(b) RESULT, CalculiX S8R, es 0.25, 7.56 rpm, dev jobs 11618743/44/47 (oml centred) and
  11618760/61/68 (oml OFFSET); `logs/check/tipmodes/t6-ccx-{oml,offset}/result.json`.
  Prestress check L8(b): PASS, all windows within 1.2 %. Local N_zz field per element vs ours (same mesh):
  corr 1.000 (0.996 at z/L 0.94-0.97), compressive means within 1-5 % -> the prestress, including the
  local TE compression, is the same in both codes (N_xx and N_xy not compared).
  Buckling, lowest positive: CCX oml 2.290, 2.572, 4.18; CCX offset 2.257, 2.547, 3.79; peak in
  09_46/09_47_LP_TE_PANEL at z/L 0.96. Ours: 0.518 (oml) / 0.541 (mid), a factor ~4.4 lower.
  Rotating frequencies, CCX vs ours: oml 0.5742/0.7047/1.636/2.094 vs 0.5742/0.7044/1.652*/2.095;
  offset 0.5449/0.6990/1.520/2.072 vs mid 0.5474/0.7007/1.537/2.074. Beam modes agree within 0.5 %
  (2F within 1.1 %; *our oml 2F is mixed with a local, beam fraction 0.828). In CCX the first local mode
  is at 3.61 (oml) / 3.60 Hz (offset); in ours, locals sit at 1.53/2.36/2.85/3.45 Hz (oml) and 2.17/2.92/2.96/3.61 Hz (mid).
  VERDICT under the L8(c) rule: UNRESOLVED. lam_ccx = 2.29 lies in (1.5, 3), and one CCX local mode at
  3.61 Hz is below 3.8 Hz. Substance: an independent shell code with the same mesh, layup and local
  prestress keeps the tip TE panels stable at rated (lam ~ 2.3), with no local mode below 3.6 Hz. Ours
  calls them buckled (0.52). The discrepancy is isolated to the element level (local K_G and/or local
  bending K of MITC4 vs S8R/C3D20R), not to the prestress, the input or the reference surface: oml vs
  offset/mid moves lam by under 5 % in both codes, which also weakens B*kappa. The global stiffening (beam
  modes) agrees. Discriminating run: ours at es 0.15 (dev 11618782). T4 already showed the lowest local
  frequency rising +34 % from 0.25 to 0.15, i.e. towards CCX.
- L14 (2026-10-19) Maintainer hypothesis: "en la punta hay mas elementos mitc3 que mitc4 es algo que deberia
  ser evaluando". REFUTED at es 0.25 (prestress npz, n_nodes per element). The mesh has 92 triangles out of
  33813 elements. There are 0 triangles in z/L 0.90-0.97, where the buckling and local modes peak; 2 of the
  509 LP_TE_PANEL elements at z/L >= 0.86; and 10 (1.3 %) at z/L >= 0.97, the tip cap. Triangles carry
  <= 0.5 % of the negative K_G energy of every local mode (oml and mid). The element-level discrepancy
  (L13) is in MITC4, not MITC3.
- L15 (2026-10-19) Ours, oml es 0.15 (dev 11618782, 2m59s, `buckle-oml-015-dev-11618782.log`): lam = 0.480,
  0.520, 0.567, 0.634, 0.706, ... with the same LP_TE_PANEL 09_42..09_47 modes at z/L 0.86-0.97 and residuals
  ~1e-10. lam_min goes 0.518 -> 0.480 from es 0.25 to 0.15: it does NOT move towards CCX (2.29). This rules
  out under-resolution of the MITC4 mesh; the ~4.4-5x gap is a formulation difference (local K_G of MITC4
  and/or local panel bending/transverse-shear K for these sandwich laminates). Gap found: no test validates
  MITC4 out-of-plane plate buckling against a closed form (grep of tests/ for buckl* or K_sigma finds only
  units/scaling tests in `tests/validation/element/test_stress_stiffened_solver.py`). Proposed T7: an
  analytic plate-buckling anchor. (i) Isotropic SS square plate, uniaxial: N_cr = 4 pi^2 D / b^2.
  (ii) Orthotropic CLT closed form. (iii) A sandwich with a soft core, with the transverse-shear correction.
  (iv) Our parked panel stiffness (K without K_G) on the TE panel vs CCX, to split K from K_G.

## Handoff 2026-10-19 (2) - T1-T4, T6 done; T7 next in a new session
Supersedes the handoff below. State: K_G alone causes the tip local modes. The prestress is identical to
CalculiX's, and the beam modes agree within ~1 %. MITC4 buckles the LP_TE_PANEL 09_42..47 panels 4-5x
earlier than CCX (lam 0.52 / 0.48 at es 0.25 / 0.15 vs 2.29). This is not the mesh, the MITC3 elements,
the input or the reference surface (L12-L15). Next session, T7 (S7):
1. Resume: this section plus L13-L15. `mem_search "mitc4-panel-buckling-gap"`. Load `.pi/skills/fem-shell/SKILL.md`.
2. Analytic anchors, each with its reference checked twice before comparison (rule L20 of
   static-dtube-before-blade.md). (i) Isotropic SS square plate under uniaxial compression,
   N_cr = 4 pi^2 D/b^2, with a mesh sequence. (ii) Orthotropic CLT closed form. (iii) A soft-core sandwich
   with the transverse-shear correction. Use `assemble_geometric_stiffness(stress_field=...)` with a uniform
   sigma (Pa) and the Arnoldi solve from `diagnose_rotating_prestress.py buckle`; scipy eigsh fails here (L11).
3. (iv) Split K from K_G on the real panel: parked panel stiffness ours vs CCX
   (`tools/t6_ccx_centrifugal_buckling.py` decks).
4. Then T5. Runs under 20 min go to sequana_cpu_dev, one at a time (Engram `config/slurm-partition-choice`).
Uncommitted, on main: `tools/diagnose_rotating_prestress.py`, `tools/t6_ccx_centrifugal_buckling.py`, this
document. Branch before committing. `scripts/known-red-full.txt` is untracked and is not this feature's.

## Handoff 2026-10-19 - T3/T4 done, T1/T2 resubmitted
Supersedes the 2026-10-18 handoff. T3 (K_G is the cause) and T4 (rated locals move +34 % with mesh,
beam modes converge) are recorded in L4/L6; the tool fix is in L5. Next session:
1. `sacct -j 11618652,11618653 -X --format=JobID,Elapsed,State`; grep
   `logs/check/tipmodes/prestress-*-1161865[23].log` for `prestress solve|wrote|Traceback|SystemExit`.
   The fixed tool has never run end to end: a new traceback there is a tool defect, not physics.
2. T1/T2: `scripts/aeroenv.sh python tools/diagnose_rotating_prestress.py analyse logs/check/tipmodes/prestress-oml-025.npz`
   (then mid). Sanity first: `rq` == `w2` per mode, `T_fe/T_ref` ~ 1 inboard. Only then: top K_G
   elements of the tip locals (z/L 0.94-0.97), Nmin < 0?, set (skin/web), h, area (slivers at the 2
   interpolated tip stations?).
3. Optional, closes T4 as a trend: oml es 0.125 rated modal (`diagnose_mode_beam_projection.py solve`, sbatch).
4. T5: record finding (closures doc, S-4 hazard: local near 2nd edge at es 0.15, beam frac 0.659).
Uncommitted: `tools/diagnose_rotating_prestress.py` (L3 tools + L5 fix), this document.
