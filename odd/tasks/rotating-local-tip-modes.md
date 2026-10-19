# rotating-local-tip-modes

Why do low-beam-fraction (local) modes appear at 7.56 rpm but not parked, with the residual at
z/L 0.94-0.97? IEA-15 official blade, es 0.25, oml and mid. Follows `odd/tasks/s4-mode-classifier.md`.
Branch: `feat/rotating-local-tip-modes` from `29897f1`; delivery: `single-pr`; checks: aeroenv + focused tests + quick.

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
- T7 [S7] analytic plate-buckling anchor for MITC4 K_G (L15 plan i-iv) - (i) DONE (L16, 3f0169d), (iii) measured with 2 strict xfails (L17, 968d20d); (iv) measured, diagnostic only (L20, SLURM 11618812), delegated worker (context-backstop trigger); (ii) next, delegated worker at context backstop - commit: (iv) this parked-panel work unit
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
- L16 (2026-10-19) T7(i) RESULT, isotropic SS plate, uniaxial compression, branch `feat/rotating-local-tip-modes`,
  `tests/validation/element/test_plate_buckling_analytic.py`, log `logs/check/tipmodes/t7-i-plate.log`.
  Reference N_cr = k pi^2 D/b^2, k = min_m (m b/a + a/(m b))^2, checked twice before the FE result was read (L20 rule):
  (A) Timoshenko & Gere k table, a/b 0.5/0.6/0.8/1/1.5/2 -> 6.25/5.14/4.20/4.00/4.34/4.00 with m 1/1/1/1/2/2: PASS;
  (B) independent FD solve of D lap^2 w = N(-w_xx) with w = lap w = 0: a/b 1 err -0.321/-0.080/-0.020 %, a/b 1.5
  -0.405/-0.100/-0.025 % at n 15/31/63 (clean O(h^2)): PASS. Tolerances fixed before measuring: monotone error over
  8/16/32, observed order in [1.5, 2.5], Richardson within 0.3 %, 32x32 within 1 %. MITC4, steel, b/h 100 (Mindlin
  correction 0.06 %), hard SS, prescribed sigma via public `assemble_geometric_stiffness(stress_field=...)`, Arnoldi
  shift-invert (residuals <= 1e-10): square N_cr err +9.25/+2.14/+0.49/+0.088 % at n 4/8/16/32, order 2.025,
  Richardson -0.045 %: PASS. a/b 1.5 at 48x32: m = 2 picked, +0.12 %: PASS. Dense eigh = Arnoldi to 1e-9 at 8x8;
  sigma_xx vs sigma_yy load identical to 1e-8. VERDICT: the MITC4 K_G formulation (N = sigma h, translational
  B_geo) is correct for a thin isotropic plate; it converges from the stiff side, so it cannot by itself make the
  blade panels 4-5x softer. With L13 (N field = CCX's), the gap moves to the laminate K: the panel bending D
  and/or the transverse-shear stiffness of the soft-core sandwich, or the K_G-K interaction at small b/h.
  Next: (iii) sandwich with shear correction and (iv) parked TE-panel K vs CCX are now the discriminating steps;
  (ii) orthotropic CLT is the cheaper intermediate.

- L17 (2026-10-19) User choice after L16: "(iii) sándwich, después (iv) (recomendado)". T7(iii) RESULT, log
  `logs/check/tipmodes/t7-iii-sandwich.log`. Section: faces 2 x 2 mm E 20 GPa, core 30 mm E 60 MPa, nu 0.3 (isotropic
  FSDT plate). Reference N_cr = min_m D K^4/kx^2/(1 + D K^2/S), checked twice: (A) numerical 3x3 Navier FSDT Hessian =
  closed form to 1e-12 for m 1-4; (B) S -> inf gives Kirchhoff, b -> inf gives Allen's strut exactly, and S from two
  definitions (own energy-equivalent quadrature 7.8895e5 vs Allen G_c d^2/t_c 7.8769e5 N/m, +0.16 %); aeroelast's
  Python `Laminate` D11 and Cs equal the own values (1e-6 / 1e-4). Choice of b: 0.5 m rejected (shear crimping, N
  decreases with m towards S, no finite critical m); b = 1.5 (m = 1, shear lowers N 25 %) and b = 0.75 (m = 2,
  shear lowers N 2.25x). MITC4: b = 1.5 N 4.041/3.976/3.960e5 at 8/16/32, order 2.025, Richardson 3.9551e5 vs FSDT
  3.1783e5 (+24.4 %): FAIL; b = 0.75 at 48x48 1.5536e6 with m = 1 vs FSDT 7.0621e5 with m = 2 (+120 %): FAIL.
  Root cause, measured: the same FSDT with S = A55 = sum G_k t_k = 3.1462e7 (39.9 x S) gives 3.95496e5 (MITC4
  +0.002 %) and 1.55271e6 with m = 1 (+0.06 %). The element is a correct FSDT plate; its composite shear input is
  wrong for sandwiches: `crates/aeroelast-py/src/assembler.rs:300-311` feeds `shear_stiffness_uncorrected()` (A55)
  to the shear block since 8cbfc0b (2026-09-28, "ADR-1 amended", Ko-Lee-Bathe 2017 p. 410 "no numerical factor").
  For a homogeneous ply that only drops the 5/6; for a soft core it replaces the core-governed S = G_c d^2/t_c by
  the face-governed sum, so sandwich panels lose their shear flexibility (unconservative local buckling and local
  modes). The blade uses this path (`build_rust_properties` -> `RustLaminate` -> `MaterialSpec::Composite`).
  Blade scope (`tools/diagnose_section_shear.py`, `logs/check/tipmodes/t7-iii-blade-shear.log`, 27 s): 676
  multi-ply sections, A55/Cs median 1.12, max 30.6, 230 above 2 (foam-core sandwiches).
  NOT the tip-mode cause: (a) direction, a too-stiff shear raises lam, and ours is 4.4x BELOW CCX; (b) the tip
  panels 09_42..09_49_LP_TE_PANEL have no core, they are 2.7-3.4 mm gelcoat/glass_triax/glass_triax, b/h ~ 100+,
  where transverse shear barely matters. Open question, not investigated: for those tip panels A55/Cs = 0.3-0.4,
  i.e. the energy-equivalent Cs is 2.5-3x ABOVE the plain integral (a near-homogeneous stack should give ~0.83).
  Tests: the two sandwich-vs-FSDT comparisons are `xfail(strict=True)` (band untouched); a characterization test
  pins MITC4 = FSDT(S = A55) and must be deleted with the fix. Fixing the section shear is a contract change (ADR-1,
  every composite result, S-2/S-4 anchors may move): needs a user decision.
  Consequence for the tip gap: with K_G (i) and the prestress N_zz (L13) cleared, and shear irrelevant for these
  thin panels, the remaining candidates are (1) N_xx and N_xy of the prestress, never compared with CCX (L13);
  (2) the bending D of the thin orthotropic triax laminate, including the per-element ply-angle correction
  (`assembler.rs` angle + angle_deg) = T7(ii); (3) curvature of the real panel (flat-plate anchors do not test it);
  (4) T7(iv) parked panel K vs CCX.

- L18 (2026-10-19) Continuation: user "continua" (twice), following the explicit order
  "(iii) sándwich, después (iv) (recomendado)". Proceed with T7(iv), parked TE-panel stiffness vs CCX;
  do not change the section-shear formulation or the physics tolerances. Route: bounded worker, because the
  runtime signalled the parent context backstop. Edit surfaces: a new panel-parity diagnostic and its focused
  tests; existing production assemblers, CCX writers and T6 tools are read-only.
  Evidence qualification: L16 validates the flat isotropic plate only; it does NOT clear K_G for a curved,
  orthotropic blade panel. L13 compared N_zz, not the full stress tensor. Thus "the gap must be in laminate K"
  was too strong: both K and the untested stress/K_G components remain hypotheses until separately measured.
  Verification plan: two reference checks before reading a panel comparison, identical extracted geometry,
  material orientation, loads and constrained boundary; report load work and residuals, not a single node.
  Physics comparisons are diagnostic unless an independently justified tolerance is available.
  Implementation: tools/t7_panel_stiffness_parity.py + tests/test_t7_panel_stiffness_parity.py (bounded worker).
  Meaningful RED: five behavioral failures after module scaffolding, then GREEN: 5 passed (2.40 s); parent
  spot check: 5 passed (2.48 s), exit 0. Ruff exit 0; quick exit 0 (75 passed, 1 xfailed, 15.92 s).
  CCX reference, before panel comparison: homogeneous axial composite strip, analytical work
  F^2 L/(2 E b h) = 2.380952381e-8 J; errors +0.66261 % and +0.31060 % on 4x16 and 8x32 meshes
  (existing 2 % membrane band). Load doubling gives 3.999998957 and 4.000003140 work ratios.
  This validates membrane/load/BC plumbing only, not curved orthotropic bending accuracy. Native assessment
  unavailable: eligible untracked paths not declared (includes unrelated known-red-full.txt); no tier or
  review approval is claimed. Returned plan retains writer self-verification, no separate verifier.
  Delivery: the feature starts at 29897f1 (local main at branch creation), not origin/main; inherited main
  history is not part of this feature. Existing slices: c6f0573 (T6 carry-over), 3f0169d + f961694 (T7i),
  968d20d (T7iii). T7(iv) is a separate local work-unit slice. Publishing/PR shape remains undecided.

- L19 (2026-10-19) First T7(iv) panel run: SLURM 11618805, FAILED exit 1 after 2m31s,
  MaxRSS 350448K, logs/check/tipmodes/t7-iv/panel-dev-11618805.log. CCX reference checks passed; actual
  panel extraction and both static solves reached the artifact-writing step (panel.npz exists). Failure:
  `TypeError: Object of type ndarray is not JSON serializable` at json.dumps(result): real NuMAD section
  metadata contains arrays not covered by the synthetic tests. No result.json was delivered; no parity
  verdict is claimed. One bounded worker correction: regression for realistic section metadata and JSON
  serialization, preserving geometry/materials/loads and already produced artifacts. Not a solver timeout.
  Correction observed: array/scalar-only JSON conversion; unsupported object types still fail visibly.
  Regression through panel(): RED 2 failed, 5 passed, exit 1 -> GREEN 7 passed (2.57 s), exit 0;
  parent spot check 7 passed (2.49 s), exit 0. Ruff exit 0; quick exit 0 (75 passed, 1 xfailed, 13.44 s).
  Original panel.npz read-only check: 100 nodes, 76 quads, 5 fixed corners; all finite, both fixed
  translations zero (CCX ~8e-25 m); load-work ratio MITC4/CCX = 0.8437892767. Preliminary artifact check
  only; a fresh end-to-end SLURM run is required to verify that the complete report is now delivered.

- L20 (2026-10-19) T7(iv) parked panel result: corrected SLURM 11618812 COMPLETED exit 0, elapsed 36s,
  MaxRSS 365652K. Evidence: logs/check/tipmodes/t7-iv/panel-dev-11618812/{result.json,panel.npz,panel.inp,
  panel.stdout}; submission: logs/check/tipmodes/t7-iv/panel-dev.sbatch (sequana_cpu_dev, 1 process, 8 GB,
  20 min cap). Revision 968d20d plus the new diagnostic; source SHA-256 is captured in the SLURM log.
  Same geometry/connectivity, centered layup, shared orientation bucketing, inboard-cut clamp including CCX
  midsides, and serialized corner-node force vector. Patch: 09_46_LP_TE_PANEL, 100 nodes, 76 quads, 5 fixed
  corner nodes; es 0.25, official IEA-15, OML geometry with OFFSET=0. Actual layup is Gelcoat 1 mm,
  glass_triax 1.0003966 mm, glass_triax 0.7109679 mm (not symmetric, no core).
  Both reference checks passed again before the panel: axial work errors +0.66261/+0.31060 % on 4x16/8x32,
  load-doubling work ratios 3.999998957/4.000003140. These validate axial reference plumbing, not curved
  orthotropic bending accuracy. The same discrete ~1 N resultant was used by both solvers; serialization
  changes the requested vector by 1.88e-7 relative (not a load mismatch between solvers).
  W_MITC4 = 0.00829268508357279 J, W_CCX = 0.009827910015117946 J, ratio 0.8437892767451503
  (-15.6211 % compliance; reciprocal directional effective stiffness ratio ~1.18513).
  MITC4 free residual 5.58e-8, energy identity error 1.10e-10, force-balance error 6.11e-9; all below
  their pre-existing diagnostic limits. No physical parity band is defined: verdict DIAGNOSTIC ONLY.
  Interpretation: this parked patch/load is stiffer with MITC4, so it does not reproduce a 4-5x softer
  bending response. It does NOT clear the whole K matrix, K_G, the full stress tensor, or the blade
  buckling mode: the isolated clamp and load do not match that eigenmode. Still open: orthotropic anchor
  T7(ii), N_xx/N_xy parity, and curved-panel stress/geometric-stiffness effects. Do not change anchors.
  Checks: final diagnostic regression 7 passed (parent 2.49 s), Ruff clean, quick 75 passed/1 xfailed;
  existing sandwich strict xfails retained, no full-suite/Rust build (no production or Rust edits).
  Native assessment unavailable (undeclared untracked paths); no native review approval is claimed.
  Risk: none for this isolated diagnostic; existing shared helpers were only read, not changed.
  Delivery workload: this coherent diagnostic is 498 new source/test lines (307 + 191), plus logbook;
  do not split code from tests or compress it. Before its commit, resolve feature delivery shape once.

- L21 (2026-10-19) User selected "Un único PR — menos recomendado" (`single-pr`). Whole future PR scope:
  rotating-local-tip-modes T1-T7 tooling, analytic anchors, evidence and closure, based on local main
  29897f1 (exclude inherited origin/main history). No chain strategy, tracker, or chain diagrams.
  Keep coherent local work-unit commits; no publishing, PR creation, merge, or review consent granted.
  The ~498-line T7(iv) diagnostic plus logbook is one unit because its shared reference/deck path and tests
  must remain together; no code compression or test-only split. After recording it, proceed with T7(ii),
  the remaining orthotropic CLT anchor, without a production section-shear change.

- L22 (2026-10-19) User correction verbatim: "pero no me has hablado nada sobre la identificacion de los problemas fisicos, yo no hago nada con que los tests pasen eso no es lo que estoy resolviendo".
  Prioritize identifying the cause of the physical instability discrepancy, not test/commit status.
- L23 (2026-10-19) Full prestress tensor, existing data only, common tangent axes (span projection and
  normal cross span). This corrects the earlier overstatement that the prestress is the same. Geometry
  mapping max centroid error 0; reconstruction reproduces ours N_zz relative L2 1.21e-16 and saved CCX
  N_zz 3.76e-8 (rounded exported ply thicknesses). CCX input is GLOBAL=YES, layer-integrated first static
  block; ours is local Pa, mapped by element frame and h. Area-weighted transverse resultants ours/CCX:
  09_42 -366.53/-45.15; 09_43 -446.35/-38.10; 09_44 -457.91/-47.10; 09_45 -412.05/-41.61;
  09_46 -430.12/-106.86; 09_47 -418.63/+12.12 N/m. Thus panel 47 changes from small mean transverse
  tension in CCX to strong compression in ours. For 46, span/shear means ours vs CCX are 753.52/928.36
  and 418.12/423.94 N/m; for 47, 421.81/721.97 and 623.05/582.20. This is a destabilizing difference,
  not just the parked stiffness difference. It is not yet proof of the 4.4x buckling-factor cause.
  Source qualification: our recovered membrane stress is A*epsilon/h, not (A*epsilon+B*kappa)/h
  (mitc4.rs:3041-3083, materials/composite.rs:53-57); the same field feeds K_G
  (assembly/assembler.rs:742-766). The actual unsymmetric tip stack makes missing B*kappa a candidate,
  but centroid-vs-Gauss averaging, frame conventions and other prestress differences remain alternatives.
  Causal experiment: same K and clamped free DOFs, replace only K_G stress with CCX, then replace only
  the transverse component (globally and within LP_TE_PANEL 42..49). Reproduce the stored baseline first,
  require eigenpair residuals, and compare the original-mode Rayleigh energies as well as new modes.
  Script: logs/check/tipmodes/t7-iv/stress-swap-dev.sbatch; one sequana_cpu_dev job, 20 min, 16 GB.
  No production code, physical tolerance, or shear correction change. First experiment SLURM 11618824:
  FAILED exit 1 after 1m58s, MaxRSS 2189940K. Baseline reproduced lambda 0.5176621298031965,
  original-mode Rayleigh 0.5176621298045456, residual 1.62e-10. CCX-stress branch returned no positive
  factors among the largest-magnitude inverse eigenvalues; this is NOT evidence of physical stability.
  Bounded numerical correction: target largest REAL inverse eigenvalues (nu=1/lambda), not magnitude,
  to seek the smallest positive factor; log spectrum and persist fixed-mode energy before extraction.
  Matrices, stresses and free DOFs unchanged. Corrected experiment pending.

## Handoff 2026-10-19 (4) - T7(iv) measured; packaging resolved in L21
Supersedes (3). Branch feat/rotating-local-tip-modes, base local main 29897f1, HEAD 968d20d.
New files tools/t7_panel_stiffness_parity.py and tests/test_t7_panel_stiffness_parity.py are verified but
not committed; the task document is modified. scripts/known-red-full.txt is unrelated and untouched.
Read L18-L20 for evidence, not the historical log. Successful result: SLURM 11618812, parked work ratio
MITC4/CCX = 0.843789; no physical PASS/FAIL band. Failed run 11618805 is preserved (JSON metadata defect,
fixed with RED/GREEN). T7 is not complete: (ii) remains open and (iii) retains two strict xfails.
Next: choose delivery shape for the oversized feature, record/commit the coherent T7(iv) unit, then
continue the orthotropic anchor without changing production section shear or widening any band.
Engram has active-scope progress/evidence records; the full historical document mirror remains pending.

## Handoff 2026-10-19 (3) - T7(i) and (iii) done; decision pending on the section-shear defect
Supersedes (2). Branch `feat/rotating-local-tip-modes` (from main; commits c6f0573, 3f0169d, f961694 and this one).
State: MITC4 K_G is correct for a thin isotropic plate (L16). The composite shear block consumes A55 instead of the
energy-equivalent Cs, which is wrong for sandwiches (L17; 230 of the 676 blade multi-ply sections), but it is NOT the
cause of the tip-panel gap (wrong direction; the tip panels are thin monolithic glass).
1. Resume: this section plus L16-L17. `mem_search "mitc4-panel-buckling-gap"`. Load `.pi/skills/fem-shell/SKILL.md`.
2. User decision pending: fix the section-shear input (ADR-1 amended, 8cbfc0b) now, as its own feature with anchor
   re-measurement (S-2, S-4), or defer it and keep the strict xfails.
3. Tip gap, cheapest first: (1) compare N_xx and N_xy per element, ours vs CCX, at z/L 0.94-0.97 (T6 data, L13 compared only
   N_zz); (2) T7(ii) orthotropic CLT plate with a glass_triax layup and a rotated ply angle (tests the per-element
   angle correction); (3) T7(iv) parked TE-panel K vs CCX. Then T5.

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
