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
- S8. "ante la duda, en orden de prioridades el ancla seria: 1- resultado analitico 2- resultado publicado en revista cientifia 3- resultado de un arbitro numerico externo" ... "No somos desarrolladores de Calculix, CCX es solo uno de los arbitros externos, que va en 3er ugar en la escala de prioridades de validacion."
- S9. "ok, vamos con tu recomendacion." (context L36: lambda1 vs CCX is a diagnostic, not an anchor; next, two analytic anchors: orthotropic plate buckling closed form, centrifugal prestress closed form through compute_prestress_field).
- S5. "El rojo de S-4 en el 1er modo de borde (5.9–6.1 % contra un límite de 5 %) es un hallazgo de precisión. No se toca la banda."

## Tasks

- T1 [S1] locate local-mode residual (elements, section, skin/web, tip refinement stations) - inline - DONE (L7) - commit: (tools, this handoff)
- T2 [S2] centrifugal prestress at 7.56 rpm: per-element membrane resultants N vs the section tension T(r) - inline + sbatch - DONE (L7) - commit: -
- T3 [S3] rated modal with K_G off, and with K_SP off - inline + sbatch - DONE (L4) - commit: -
- T4 [S4] es 0.25 vs 0.15 local-mode frequencies and beam fractions - sbatch - DONE (L6) - commit: -
- T6 [S6] (CLOSED as DIAGNOSTIC by S9, L36: lambda1 3.224 vs CCX 2.29 documented and localised in TE_FLAT slivers, L35; not an anchor) independent shell arbiter for the tip local modes: our linear buckling factor + CalculiX S8R CENTRIF -> *BUCKLE and *FREQUENCY PERTURBATION at 7.56 rpm, same mesh es 0.25 - worker + sbatch - commit: -
- T7 [S7] analytic plate-buckling anchor for MITC4 K_G (L15 plan i-iv) - (i) DONE (L16, 3f0169d), (iii) measured with 2 strict xfails (L17, 968d20d); (iv) measured, diagnostic only (L20, SLURM 11618812), delegated worker (context-backstop trigger); PAUSED by user (L24): analyze stress-swap 11618827 first, then (ii); (ii) orthotropic closed form DONE (L37, all predeclared bands met) - commit: (iv) af7b8fe, (ii) 3da1cb8
- T8 [S8, S9] analytic anchor for the centrifugal prestress recovery: compute_prestress_field (the path fixed in 41e9a16) vs a closed form (e.g. rotating prismatic member, sigma(r) = rho omega^2 (R^2 - r^2)/2, plus a case where B kappa is non-zero), tolerance predeclared in the docstring - inline - DONE (L38, all predeclared bounds met) - commit: (this)
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

- L24 (2026-10-19) Explicit session stop, user verbatim: "haz el cierre para analizar resultados en otra sesion q ya esta esta con el contexto muy largo". No further physical result analysis in this session.
- L25 (2026-10-19) During closure, corrected stress-swap job 11618827 COMPLETED, wrapper exit 0
  after 3m22s. Results deliberately UNREAD/UNINTERPRETED, not a physics PASS. Artifacts:
  logs/check/tipmodes/t7-iv/stress-swap-11618827/{result.json,stress-fields.npz}; stdout:
  logs/check/tipmodes/t7-iv/stress-swap-11618827.log. Reproducible submission script:
  logs/check/tipmodes/t7-iv/stress-swap-dev.sbatch (ignored local artifact, not tracked source).
  Local work-unit commit af7b8fe holds the parked-panel diagnostic, its 7 focused tests and L18-L23.
  It is a coherent 612-addition/1-deletion unit, accepted under the user's single-pr delivery choice;
  do not count inherited local-main commits as part of this feature. No push or PR created.
- L26 (2026-10-19) User request verbatim: "Analizá los resultados de stress-swap-11618827: cuánto de la brecha
  de pandeo explica el pretensado completo y cuánto la componente transversal, manteniendo K fijo. Priorizá la
  causa física; no modifiques la formulación hasta demostrarla."
  Stress swap 11618827 (same K/free DOFs; residuals <=1e-9 on low modes; first positive lambda / original-mode
  Rayleigh / original-mode destabilizing energy vs ours): ours 0.5177/0.5177/1.00 (baseline reproduced);
  ccx_all 3.189/7.17/0.072; ccx_transverse_only 2.501/6.62/0.078; ccx_transverse_tip_only (09_42..49
  LP_TE_PANEL) 0.559/0.613/0.84. Gap 0.518->2.29 is x4.42; stress alone x6.16 (overshoots), so the operator
  residual (K/K_G discretisation) is x0.72, i.e. ours stiffer, consistent with W ratio 0.844. Transverse N is
  ln-share 87 % of the stress effect (x4.83 of x6.16); span+shear x1.28. Tip-panel-only swap x1.08: the mode
  migrates (peak z/L 0.864->0.845); the spurious compression is distributed over LP_TE_PANEL and LP/HP_TE_REINF,
  z/L 0.7-1.0 (area-mean N_tt z/L>=0.8 ours/CCX: LP_TE_REINF -1108/-62, HP_TE_REINF +499/+47, LP_TE_PANEL
  -284/-35, HP_TE_PANEL +59/-5; symmetric SW -317/-326). Mismatch tracks laminate asymmetry e2=B22/A22: SW 0,
  tip TE_PANEL e2/(h/2)=0.27, tip TE_REINF 0.38 (1 mm gelcoat + 0.3-0.6 mm biax).
  CCX self-split (logs/check/tipmodes/t7-iv/bkappa_ccx.py, existing ply Gauss stresses only, 3D compliance
  ORI_P0_0, exact plane-stress Q in the common (s,t) axes, linear fit eps0+z*kappa): Gauss order zeta-slowest
  (fit residual median 4 % vs 12 % alternative, closure 0.38 vs 0.79 overall); tip TE regions close to
  5-10 %, fit 1-7 %. CCX A*eps0 / B*kappa at z/L>=0.8: LP_TE_PANEL -235/+220, LP_TE_REINF -641/+627,
  HP_TE_PANEL +82/-91, HP_TE_REINF +747/-739. Ours vs CCX A*eps0 area-L2 difference at z/L 0.8-1.0 is
  ~2.6x smaller than vs CCX full.
  Causal swap SLURM 11618888 (bkappa-swap-dev.sbatch, K fixed, COMPLETED 0:0, 4m48s): ccx_minus_own_Bkappa
  0.848 (Rayleigh 0.887, energy 0.58); same but tt only 0.855; ours_plus_ccx_Bkappa 0.864 (Rayleigh 1.06,
  energy 0.49). Interpretation as a 2x2 factorial {membrane field ours/CCX} x {B*kappa absent/CCX}:
  0.518, 0.864 / 0.848, 3.189. Removing CCX's own B*kappa (almost entirely its tt part) collapses lambda
  3.19->0.85 (73 % of the ln stress effect): missing B*kappa is DEMONSTRATED as a necessary cause. It is NOT
  sufficient: adding CCX's B*kappa to ours only reaches 0.864, because our A*eps membrane field also differs
  in the TE reinforcement (LP_TE_REINF -1108 vs CCX A*eps0 -641; z/L 0.7-0.8 large area-L2 difference) and
  the two effects interact strongly (each alone x1.6, together x6.2). Caveats: ccx_minus_own_Bkappa contains
  the 5-10 % split residual; NB uses CCX's kappa, not ours; cb*E2 term not quantified. No production change.

- L27 (2026-10-19) User: "continua". Our own consistent resultant, K fixed, SLURM 11618900 (COMPLETED 0:0,
  5m21s; logs/check/tipmodes/t7-iv/consistent-swap-dev.sbatch, consistent-swap-11618900/{result.json,fields.npz}).
  Prestress source reproduced exactly (sigma rel. error 0, static residual 5.2e-10). Local ABD rotation sign +1
  verified to <=1e-8 on all MITC4 quads; the 92 MITC3 triangles (0.13 % area) fail for both signs and keep
  our sigma (worker jobs 11618895 script bug, 11618898 frame gate). N = A eps_m + B kappa + D E2 per
  mitc4.rs W blocks (kappa=(2/h)e^b1, E2=(4/h^2)e^b2). Area-mean N_tt z/L>=0.8, ours -> ours_consistent vs
  CCX full: LP_TE_REINF -1108 -> -68 vs -62; LP_TE_PANEL -284 -> -37 vs -35; HP_TE_REINF +499 -> +50 vs +47;
  HP_TE_PANEL +59 -> -2.5 vs -4.9; LP_LE -73 -> -98 vs -99; LP_SPAR 10 -> 26 vs 31; SW unchanged (B=0).
  K-fixed lambda (peak z/L, residual <=1e-10): ours 0.518 (0.864); ours_A_plus_Bkappa 3.224 (0.96), cf. ccx_all
  3.189; ours_consistent 2.849 (0.79), next 3.192 (0.96); ours_consistent_tt_only 2.464 (0.961), cf.
  ccx_transverse_only 2.501. Original-mode destabilizing energy drops to 7-9 %.
  CONCLUSION: the buckling gap is caused by the K_G prestress recovery sigma = A eps_m / h, which omits the
  membrane-bending coupling B kappa (and W_02 D E2) that the stiffness itself contains; in the thin
  gelcoat-dominated unsymmetric tip TE laminates this turns ~zero transverse resultant into strong compression.
  Our own consistent resultant reproduces CCX N_tt region by region and moves lambda from 0.518 to 2.85-3.22
  (CCX 2.29). The remaining factor 0.71-0.80 is the operator (K/K_G discretisation) difference, ours stiffer,
  consistent with W_MITC4/W_CCX 0.844; separate question. CORRECTS L26 "not sufficient": the 0.864 of
  ours_plus_ccx_Bkappa mixed our A eps_m with CCX's split (ours B kappa LP_TE_REINF +1040 vs CCX +627);
  the A eps / B kappa split is model-specific, only the total N is comparable. Open details: D E2 is
  negligible except on TE_FLAT (0.03 m2, -2679 N/m) and moves the first consistent mode to z/L 0.79 (2.85 vs
  3.22) - check whether that mode sits on degenerate TE_FLAT elements; MITC3 recovery not examined.
  Still no production change; a fix touches crates/ recovery feeding K_G and will move S-4 and every
  K_G-based anchor - user decision required.

- L28 (2026-10-19) User chose "Antes ubicar el modo de z/L 0.79". Worker job 11618908 (COMPLETED, 3m34s;
  11618903 failed on a JSON uint64 script bug), locate-mode-dev.sbatch, locate-mode-11618908/{result.json,
  modes.npz}; K fixed, reused consistent-swap-11618900 fields. consistent mode 1 = 2.8485 (residual 2.4e-11),
  peak z/L 0.790 on 09_38_LP_TE_PANEL/10_38_LP_TE_REINF; signed energy shares LP_TE_REINF 0.93,
  HP_TE_REINF 0.14, TE_FLAT -0.069 (stabilising), TE_FLAT |translation|^2 5.7 %. Removing D E2 ONLY on
  TE_FLAT: mode 1 = 3.2246 at z/L 0.960 (LP_TE_PANEL 0.95), equal to dropping D E2 everywhere (3.2244).
  So the 2.85 TE_REINF mode is not on TE_FLAT but its factor is set by how much TE_FLAT restrains the
  TE_REINF. TE_FLAT is a degenerate sliver strip: median h 24.6 mm, area/element 2.8e-4 m2 (~1 mm wide for
  0.25 m long elements), so h/width ~25 and its shell resultants are not physical in either model. Area-mean
  N_tt ours/consistent/CCX at z/L 0.6-0.8: HP 13035/1903/-168, LP 5929/1073/-171 N/m (Nss agrees within
  ~7 %). Correction to L27: TE_FLAT is 9.5 m2 over the whole span; the 0.03 m2 figure is z/L>=0.8 only.
  Reading: the B kappa conclusion stands; the 2.85-3.22 bracket is the sensitivity to the TE_FLAT sliver
  strip (a meshing/modelling artefact shared with CCX in different form), not to the panel physics.

- L29 (2026-10-19) User verbatim: "podemos hacerlo aqui  asi que corrige el codigo de produccion  luego
  probamos todo en otra sesion" and "antes de ejecutrar todo hay que validar que los cambios efectivamente
  corrijen lo que estamos evaluando". Production fix (worker): K_G prestress = (A eps_m + B kappa)/h at the
  element centre, local frame; D E2 deliberately excluded (L28); output/postprocess stresses unchanged.
  mitc4.rs/mitc3.rs compute_prestress_membrane_stress, assembler.rs compute_prestress_field (used by
  assemble_geometric_k_from_disp), PyO3 compute_prestress_field, stress_stiffened.rs and
  stress_stiffened_dynamic.py K_G sources; tests/test_kg_prestress_coupling.py; docs/formulations/
  shell-elements.md. RED observed (missing API / E0425); GREEN 6 Rust + 6 Python; maturin develop --release;
  plate-buckling anchor 24 passed 2 xfailed (unchanged strict xfails); check.sh quick 75 passed 1 xfailed.
  VALIDATION of the production path, SLURM 11618946 (validate-fix-dev.sbatch, COMPLETED, 1m45s, uncommitted
  tree on 54a6386): field new vs diagnostic N_AB/h on quads 8.2e-17, SW (B=0) new vs old 1.3e-16; triangles
  change 54 % (now get B kappa; diagnostic had kept them). assemble_geometric_stiffness(omega, free_dofs):
  lambda1 3.2243868 (residual 7.5e-11, peak z/L 0.960) vs diagnostic 3.2243628 (rel 7.4e-6, the triangles);
  next 4.067, 4.300 (z/L 0.79), 5.247. The production fix reproduces the demonstrated correction.
  Flag: Rayleigh of the OLD first mode is 9.37 vs 7.29 in the diagnostic - only the triangle change differs;
  not investigated (old mode is no longer the critical one).
  Not done (next session, user's plan): full suite, S-4 rotating modal vs OpenFAST, S-5/S-6, T7 vs CCX 2.29,
  re-derive known-red. tools/diagnose_rotating_prestress.py:111-163 still builds K_G from membrane-only sigma
  (diagnostic baseline 0.5177 tool) - update it before reusing it; test_stress_stiffened_solver.py:285-301 and
  tools/debug_kg_magnitudes.py:88 use the old pattern (isotropic / print only).

- L30 (2026-10-19) User verbatim: "prepara el handoff de esta sesion para comenzar en la proxima, pero de igual
  forma creo que antes de correr la suite ocmpleta deberiamos correr la parte de la suite que nos trajo a estos
  diagnosticos, luego toda la suite". Session closed; Handoff (9) orders the targeted K_G suite before the
  full suite.

- L31 (2026-10-19) User verbatim: "mi premisa es que todos son codigos numericos, diferentes, evaluando el mismo problema fisico. Entonces si todos coinciden dentro de una banda la fisica esta dentro de lo esperado. 5% es una restriccion mia, para evaluar que tan lejos estamos pero 5-6% tampoco es mucha diferencia considerando que son metodos diferentes" and "me preocupa mas esto - Abierto: lambda1 sigue 41 % por encima de CalculiX (factor de operador 0,844)" and "vamos a ejecutar la suite completa para ver que se movio y despues escribe el handoff". Step 1 of Handoff (9) executed (results below); full suite submitted as SLURM 11619005 (sequana_cpu_bigmem, scripts/check_sbatch.sh full, revision c197257+dirty: only tools/diagnose_rotating_prestress.py differs). Rebuilt the .so first (installed copy predated 41e9a16): maturin develop --release OK.
  STEP-1 RESULTS (logs/check/kgrevalid/, baselines in logs/check/modecls and logs/check/tipmodes untouched):
  Targeted tests, small SLURM jobs: dev1 11618971 38 passed 2 xfailed (kg_prestress_coupling, plate_buckling, stress_stiffened_solver; the 2 xfails are the strict T7(iii) ones); dev2 11618985 117 passed (omega_rebuild_predicate, rotor_centrifugal_branch, t7_panel_parity); dev3 11618990 53 passed 6 skipped (rotor_rust_parity, inertial parity x2; skips = preCICE mock, tasks 6.5-6.10); s0defl 11618973 11 passed (S-0 + deloading); campbell 11618974 8 passed. S-4 11618972 (test uses skin_reference="mid", es 0.25): 2 failed 1 passed: 1st edge parked 0.6951 vs OpenFAST 0.7390 = 5.9 %, rotating 7.56 rpm 0.7006 vs 0.7461 = 6.1 %; the same numbers as before the fix, parked is K_G-independent; the test stops at the first failed assert so 2F rotating was not printed. Band NOT widened. No anchor turned red because of 41e9a16.
  Rotating modes (oml es 0.25, 7.56 rpm; modes-oml-025-kgfix.npz, SLURM 11618948) old K_G -> new: 1F 0.57423 -> 0.57420, 1E 0.70440 -> 0.70436, 2F 1.6521 (beam fraction 0.828) -> 1.6366 (0.999), 2E 2.0947 -> 2.0907, 3F 3.1232 -> 3.0840. The four low-beam-fraction local modes (1.5325, 2.3636, 2.8491, 3.4524 Hz at z/L 0.94-0.97) DISAPPEAR; first non-beam mode now 3.716 Hz (z/L 0.84, also present parked at 3.834). Parked frequencies identical to the baseline (control). vs CalculiX (L-T6, same mesh/layup): 1F 0.5742 vs 0.5742, 1E 0.7047 vs 0.7044, 2F 1.636 vs 1.6366 (+0.04 %), 2E 2.094 vs 2.0907, first local 3.61 vs 3.716 (+2.9 %).
  Buckling (tools/diagnose_rotating_prestress.py now uses compute_prestress_field; solve 11618949 ok, buckle 11618981 with --sigma 3.0): lambda 3.22439 (z/L 0.960, 09_47_LP_TE_PANEL, beam fraction 0.115, residual 2.7e-11), 4.0669 (0.965), 4.2997 (0.790), 5.247 (0.866); production 3.2243868, so diagnostic = production. CCX 2.290/2.572/4.18. Buckle failures 11618949 (sigma 0.5) and 11618958 (sigma 0.3): the shift must be above lambda1 now, otherwise ARPACK returns only the negative (tension) branch (-2.1..-0.53). Outputs: prestress-oml-025-kgfix.npz, buckle-oml-025-kgfix.npz.
  READING: the cause (K_G without B kappa) is fixed and demonstrated against two independent codes; beam modes agree with CCX within 0.2 %; 1E vs ElastoDyn stays 5.6-6.1 % but ours equals CCX and BeamDyn parked (0.6986 vs 0.6946/0.6995), so the residual is ElastoDyn-side (ElastoDyn parked 1E 0.7403 vs BeamDyn 0.6946). User: 5-6 % between different codes is acceptable context, 5 % stays his own yardstick, not widened. Main open concern (user): lambda1 3.22 vs CCX 2.29, +41 %.

- L32 (2026-10-19) User: fix all infrastructure problems; TMPDIR under $SCRATCH/tmp (rule to AGENTS.md at close); relaunch all, dev where possible.
  Full-suite SLURM 11619021 was NOT a valid measurement: 152 pytest setup ERRORs (the inherited per-login-session TMPDIR
  $SCRATCH/tmp/gpZFUcEV vanished mid-run) that check.sh did not count, plus 3 spurious NEW-FAIL. Fix e4473ee: check.sh owns
  TMPDIR=$SCRATCH/tmp/check-<mode>-<jobid>, pytest --basetemp there, and any ERROR fails the gate (exit 1). Rerun in dev (11619071,
  4 tasks): test_blade_mass_mapping x2 and test_blade_parity_ccx[tip_axial_compression] PASS (they were tmp-induced);
  test_blade_mesh_generation x2 stay red (9357 vs 9277 nodes) and are in the known-red registry recorded at 52e9bba, i.e. red
  before 41e9a16. Valid reading of 11619021 for the rest: S-1/S-2/S-4/S-6/V-02 reds identical to the registry (S-4 1st edge 5.9 %/6.1 %
  vs 5 %); test_blade_modal_gap_converges_with_mesh now passes (stale entry). No --record done. Fresh full suite: SLURM 11619070
  (sequana_cpu, --ntasks=4, snapshot logs/check/git-snapshot-20261019T173239, revision e4473ee clean) - READ ITS RESULT.
  Step B(i)/(ii) K-vs-K_G split on the CCX buckling mode shapes (logs/check/tipmodes/t7-v/{kenergy.py,kenergy-dev.sbatch},
  SLURM 11619079, results kenergy-11619079/{ours.json,split.json}): translations of CCX mode k prescribed in CCX (RF -> phi.K_ccx.phi) and in
  ours (rotations condensed). Calibration, CCX parked-frequency mode 1 (beam): K_ours/K_ccx = 1.009 -> plumbing valid.
  CCX buckle modes 1/2/3: K_ours/K_ccx = 2.93 / 7.74 / 4.03; G_ours/G_ccx (G_ccx = K_ccx/lambda_ccx) = 0.755 / 0.765 / 0.687;
  our Rayleigh on the CCX shape 8.89 / 26.0 / 24.5 vs lambda_ccx 2.29 / 2.57 / 4.18. READING: on the CCX local mode shapes our K is
  3-8x stiffer than CCX's, but our eigenmode relaxes to a different, cheaper shape (lambda 3.22), so the K ratio of a non-eigen shape
  does NOT decompose the 1.41 eigenvalue ratio. Not yet localized (which elements/regions carry our excess energy: TE_FLAT slivers? TE_REINF? midside/rotation
  consistency of the CCX field?). Hypotheses (iii)-(v) untouched. No production change, no band touched.

- L33 (2026-10-19) User: the 3-8x K ratio of L32 does not fit what was established before (parked panel W ratio 0.844, beam modes
  within 0.2-1 %, lambda ratio 1.41). Agreed: L32 is NOT evidence about K. Likely measurement defects, untested: (a) the CCX field is
  S8R (midside nodes carry the shape) but our MITC4 uses corner nodes only, so sampling it at corners may inject spurious bending
  energy on locally short-wavelength modes (z/L 0.96, 1-2 elements wide); (b) all three buckle modes give K_ccx ~1.01-1.02 whatever the
  shape (normalisation artefact or RF sum over a subset); (c) rotations condensed in ours vs CCX expanded-solid kinematics; (d) energies
  are not compared on the same element set. Do not quote the 2.93/7.74/4.03 as a finding. First action next session: self-consistency
  test of the procedure (prescribe OUR eigenmode and a smooth analytic field in both codes; require K_ours/K_ccx to reproduce the
  0.84-1.0 seen before) before any per-element localisation.

- L34 (2026-10-19) Step B self-test of the prescribed-field K-energy procedure (L33), BEFORE any localisation.
  Script logs/check/tipmodes/t7-v/kselftest.py (+kselftest-dev.sbatch), SLURM 11619083 (sequana_cpu_dev, --ntasks=16, 3m15s, 0:0;
  tree 9a712de, no code change; the log lacks the revision because git is not on PATH before aeroenv.sh), results
  selftest-11619083/{ours.json,split.json}. Fields: OUR1 = production buckle mode 1 (logs/check/kgrevalid/buckle-oml-025-kgfix.npz,
  lambda 3.22438679), BEAM = constant-curvature flap bending with plane sections (u_y = k s^2/2, u_z = -k y s), CCXB1-3 / CCXF1 = CCX
  modes from the frd INCLUDING midside nodes. CCX midside variants: free (corners only, midsides relax = L32), edge (midside = mean of
  its two corners), full (exact midside). Ours: corner translations prescribed, rotations condensed. Criteria fixed before running:
  P1 Rayleigh of our mode = 3.224387; P2 CCX full = 1 +-1 % (normalisation); P3 BEAM ratio 1.00 +-0.02; P4 OUR1 ratio in [0.84, 1.19]
  (L20 compliance 0.844 <-> stiffness 1.185).
  RESULTS: P1 PASS 3.2243867889 vs npz 3.2243867889 (OUR1-cond = OUR1-own to 1e-12: K_G acts on translations only, so our eigenmode
  rotations are already K-condensed). P3 PASS: K_ours/K_ccx BEAM free 1.0016, full 0.9997, edge 0.9995. P4 PASS: OUR1 free 0.9938
  (ours 8.8376 vs CCX 8.8926): on OUR mode shape both codes agree on K within 0.6 %. P2 FAIL by +3 %: CCX full = 1.0305 / 1.0346 /
  1.0329 for B1/B2/B3 (free 1.0097 / 1.0238 / 1.0227); a ~3 % error of OUR round-trip procedure (frd 6 digits, midsurface translations re-imposed as BCs in a second CCX run, RF summed) remains unexplained - NOT a CalculiX error (corrected in L36); the "=1" reference is itself an inference (1+1/lambda argument), not confirmed in CCX documentation (frd 6-digit
  rounding or 2D-node BC vs expanded-solid field; not chased: it cannot touch a factor 3-8). CCXF1: free 13.41, full 13.60, ours 13.535.
  'edge' is NOT a valid comparator on curved geometry: CCX on OUR OWN mode with chord-averaged midsides costs 113 vs 8.89 (x12.7) -
  averaging corner displacements of a rotation-like field puts the midside off its true displacement and injects membrane strain.
  Q3 answered: K_ccx(free) ~ 1.01-1.02 on the buckle modes is the normalisation of CCX buckling eigenvectors, phi^T K0 phi = 1 (ARPACK
  buckling mode, B = K0; if the base stress were in B the value would be 1 + 1/lambda = 1.44 for mode 1). Not a subset artefact: the
  corner-only RF sum equals the all-node sum (free midside |RF| <= 8.2e-7); extra_nodes 66471 = the midside count, all with zero RF.
  VERDICT: the 'free' procedure (min K energy given corner translations, interior kinematics relaxed in each code) passes the user's
  self-test (our mode 0.994, smooth field 1.002); the L32 ratios 2.93 / 7.74 / 4.03 are therefore a valid measurement of THAT quantity.
  Reading (hypothesis, not yet localised): CCX's buckle modes are cheap for S8R given corner values but 3-8x dearer for MITC4, while our
  mode is equally priced by both; CCX-edge/CCX-full = 9.9 / 8.0 / 19.9 shows the CCX modes carry sub-element (midside) content.
  Next: per-element localisation, SLURM 11619087 (klocal-dev.sbatch, stages elem/local in kselftest.py).

- L35 (2026-10-19) Step B localisation (allowed: the 'free' procedure passed its self-test, L34). SLURM 11619087
  (klocal-dev.sbatch, sequana_cpu_dev, --ntasks=6, 2m41s, 0:0); split/local re-run on the login node after a parser fix (the ELSE block
  broke the RF parser; heavy stages were fine). Output logs/check/tipmodes/t7-v/klocal-11619087/{ours-elem.npz,split.json,local.json}.
  Method: same corner-translation field in both codes ('free'); CCX per-element *EL PRINT ELSE (= 1/2 uKu, doubled), ours per element
  from the contiguous assemble_k COO blocks (asserted: block count, each block touches exactly its own nodes), matched by corner-node set.
  Closure checks: ours sum-per-element = global to 1e-12; CCX 2*sum ELSE = sum RF.u (B2 1.0238 both).
  Controls: BEAM every set 1.001-1.018 and every z/L band 0.992-1.015; OUR1 total 0.994, TE_PANEL 1.056, SPAR 1.057, TE_FLAT 1.515,
  TE_REINF 0.527, SW 0.716 (compensating moderate differences on our own mode).
  CCX buckle modes, K_ours/K_ccx by set (share of ours / share of CCX):
    B1 total 2.93: TE_FLAT 116 (0.395 / 0.010), TE_PANEL 1.67 (0.485 / 0.855), TE_REINF 2.24; z/L 0.94-0.97 carries 99 %.
    B2 total 7.74: TE_FLAT 334 (0.487 / 0.011), TE_PANEL 3.59 (0.404 / 0.871), TE_REINF 5.60, SPAR 8.6, LE 1901 (1 % of ours).
    B3 total 4.03: TE_FLAT 23 (0.579 / 0.102), TE_REINF 2.27 (0.298 / 0.531), TE_PANEL 1.33; z/L 0.8-0.9 carries 99 %.
  Top-excess elements are all LP_TE_FLAT (B1: 947/973/921 at z/L 0.962; B2: 1025/1055/999 at 0.959; B3: 4271/4303 at 0.838).
  Geometry: those TE_FLAT elements are slivers 1.25 x 126 mm (B3 region 1.68 x 141 mm), laminate thickness 1.0 mm (width ~ thickness,
  aspect ratio ~100), flat (warp 0). The CCX modes impose on them eps22 ~ 2-3e-4, in-plane shear ~ 1e-3 and an across-width normal
  relative displacement dn/W ~ 0.01-0.03 (a rotation of the strip); OUR1 imposes the same order (7.6e-4 / 5.4e-4 / 2.5e-2), yet there
  ours/CCX on TE_FLAT is only 1.5. A membrane estimate of uKu for such an element (E~2e10, h 1e-3, eps 1e-3, area 1.6e-4 m^2) is ~3e-3;
  ours carries 0.77 on element 1025 (B2): the excess is not membrane energy of the strip.
  READING: the 3-8x on CCX shapes is localised in the TE_FLAT slivers (40-58 % of our energy vs 1-10 % of CCX's), plus a secondary 1.3-3.6x
  in TE_PANEL/TE_REINF next to them. Which part of our element (MITC4+D drill-membrane of Ko, Bathe & Zhang 2025, transverse shear, or
  bending) carries it is NOT determined; nor is which code is right on a 1.25 mm wide, 1 mm thick shell (neither is a valid shell there).
  This is hypothesis (iv) TE_FLAT (lambda1 bracket 2.85-3.22, L29-L30), now with a localised target. No production change, no band touched.

- L36 (2026-10-19) User: anchor hierarchy S8 (analytic > published > external numerical arbiter; CCX is only a level-3 arbiter) and
  "no entiendo como puede haber una diferencia de un 3% en un codigo externo que no tocamos": my L34 wording "CCX-side systematic" was
  wrong; the 3 % is the uncertainty of our frd -> BC -> RF round trip, corrected in L34. OpenSees asked and re-assessed: not viable for the
  blade (L8: no +-45 deg plies in this build, singularity for GLIBC), level 3 like CCX, ShellMITC4 shares our element family (ASDShellQ4
  = ShellMITC4 to 4 digits in the probe), no native centrifugal load, and (my understanding, unverified here) no linear K_G eigen analysis;
  above all the discrepancy lives in 1.25 mm x 1 mm slivers where no shell element is valid. Decision S9: lambda1 vs CCX is a DIAGNOSTIC
  (T6 closed as such); the roadmap keeps the qualitative result (lambda1 > 1 in both codes: no local buckling at rated, the pre-fix 0.518
  was the spurious one). Next: level-1 anchors for the ingredients - T7(ii) orthotropic plate buckling closed form and T8 centrifugal
  prestress recovery closed form. Existing level-1 coverage: tests/validation/element/test_plate_buckling_analytic.py (isotropic flat
  SS plate, prescribed stress, Timoshenko & Gere); tests/validation/rotor/test_rotor_inertial.py checks nodal m omega^2 r, NOT the
  recovered prestress.

- L37 (2026-10-19) User verbatim: "leé el Handoff (13) de odd/tasks/rotating-local-tip-modes.md. Empezá por el log de 11619070 contra
  scripts/known-red-full.txt y después seguí con la T7(ii). Si el contexto crece, dejá la T8 para otra sesión."
  Step A: SLURM 11619070 (sequana_cpu, e4473ee) still RUNNING at 1h05 and only 25 % (logs/check/full-20261019T173319.log; the
  slurm log is empty until exit), limit 4h00: no comparison with known-red-full.txt yet. Failures visible so far: 2 near 5 %, ~14 near
  25 % (names unknown until the summary). Risk: the run may hit the 4 h limit.
  T7(ii) DONE: tests/validation/element/test_plate_buckling_analytic.py, specially orthotropic SS cross-ply [0/90/90/0], carbon/epoxy
  (E1 140, E2 10, G12 5 GPa, nu12 0.3), h 4 mm, b/h 250, prescribed sigma. Reference N_cr = pi^2/b^2 [D11 (mb/a)^2 + 2(D12+2D66) +
  D22 (a/(mb))^2] checked twice before the FE: A (isotropic reduction = T&G to 1e-12, m switch at (D11/D22)^(1/4) sqrt(m(m+1)), wide
  column, Navier FSDT S -> inf), B (FD of the orthotropic PDE, 0.2 % at n=31, O(dx^2)). Section: aeroelast Laminate D = own CLT, B = 0;
  shear lowers N_cr < 0.1 % for min(Cs, A55-uncorrected) (asserted). FE (MITC4) vs CLT: x load (m=1) 8/16/32 = +2.47/+0.58/+0.12 %,
  observed order 2.03, Richardson -0.033 % (band 0.3 %); y load (m=2 picked, lam2/lam1 1.086) 48x48 +0.15 % (band 1 %);
  [90/0/0/90] under y == [0/90/90/0] under x to 1e-8 (ply angle follows element e1). Reading: MITC4 K_G + composite D is a correct
  Kirchhoff orthotropic plate where the shell is thin; no band touched. Two reference-check tolerances were mine and numerical
  (S = 1e20, then 1e13, singular/round-off in the 3x3 Navier solve), fixed to S = 1e11, rel 1e-6; no FE band was changed.
  Checks: file 25 passed, 2 xfailed (T7(iii) strict xfails unchanged); scripts/check.sh quick OK (75 passed, 1 xfailed).


- L38 (2026-10-19) User verbatim: "vamos con estas dos tareas en paralelo / 1. Suite: "leé el Handoff (14) y hacé el paso A". Compara
  el resultado con known-red-full.txt. Si se pasó de las 4 h, hay que reenviarlo a sequana_cpu_bigmem con más --time. / 2. T8: "leé el
  Handoff (14) y hacé T8"."
  Step A CLOSED: SLURM 11619070 (sequana_cpu, e4473ee) COMPLETED 02:15:27 of 4h00, no resubmission. pytest 20 failed, 1357 passed,
  38 skipped, 15 xfailed, 3 xpassed (names not in the log), 0 ERROR lines; check.sh full OK. The 20 failures equal the 20 entries
  of scripts/known-red-full.txt (52e9bba) one-to-one: none outside the registry, none stale (test_blade_modal_gap_converges_with_mesh
  still red, so it is not stale). No 54a6386 rerun needed; no --record.
  T8 DONE: tests/validation/element/test_prestress_centrifugal_analytic.py. Flat strip r0 2 m -> R 10 m, b 0.5 m, root clamped, spin
  about y (body force in plane), nu = 0 so the uniaxial state N_x = m w^2 (R^2 - x^2)/2 is the exact plane-elasticity solution; M_x = 0
  (no transverse load) so N is set by equilibrium for any layup. Discretisation predicted before measuring: centroid-lumped load gives
  N_i = N(x_c) - m w^2 L^2/8 exactly. Measured (isotropic, symmetric, unsymmetric [soft/stiff], 20 and 40 elements): normalised error
  4.1667e-4 / 1.0417e-4 = the predicted bound, ratio 4.00; discrete identity <= 2.1e-11; sigma_yy, sigma_xy <= 1.7e-13. B kappa:
  output A eps/h over prestress = A11 D11/(A11 D11 - B11^2) = 2.008299 to 1e-6 (unsymmetric), 1 (symmetric): the pre-41e9a16
  recovery reads this tension 100 % high. Not covered: MITC3 (uneven load split across the width), nu != 0 (not uniaxial).
  Harness slip on the first run (assembler element proxies carry no nodes), fixed; no bound was touched.

## Handoff 2026-10-19 (15) - START HERE: T5 closure
Supersedes (14). Step A and T8 DONE (L38). Remaining: T5 (record the finding in docs/validation_closures.md, S-4 docstring only if
a number moved: the 11619070 S-4 reds are the registered ones), TE_FLAT sliver decision (L28), old-mode Rayleigh flag 9.37 vs 7.29
(L29). No push; branch feat/rotating-local-tip-modes.

## Handoff 2026-10-19 (14) - superseded by (15): T8, then read full suite 11619070
Supersedes (13). T7(ii) DONE (L37, 3da1cb8). Session closed by the user at 19:20 (going home); 11619070 was at 60 % at
1h47 (4 h limit, likely to finish). Resume manually once squeue shows it gone. Step A: read
logs/check/full-slurm-11619070.log tail vs scripts/known-red-full.txt (untracked); if it timed out, resubmit with a longer --time on
sequana_cpu_bigmem; outside-registry failures need a 54a6386 worktree rerun of that test alone; no --record unless justified.
Then T8: compute_prestress_field vs a closed-form centrifugal stress (rotating prismatic member sigma(r) = rho omega^2 (R^2 - r^2)/2,
plus a case with non-zero B kappa: offset reference surface or unsymmetric laminate), tolerance predeclared in the docstring.

## Handoff 2026-10-19 (13) - START HERE: level-1 anchors T7(ii) and T8 (decision S9)
Supersedes (12). Read S8-S9, L36 and the T7/T8 lines. Step A unchanged: full suite SLURM 11619070 was still RUNNING at 1h01 at close;
read logs/check/full-slurm-11619070.log tail vs scripts/known-red-full.txt (untracked); outside-registry failures need a 54a6386
worktree rerun of that test alone; no --record unless justified. Then T7(ii): orthotropic SS plate under uniaxial compression,
N_cr = pi^2/b^2 [D11 (m b/a)^2 + 2 (D12 + 2 D66) + D22 (a/(m b))^2] (Whitney, check the form and its reference twice before reading
the FE, as test_plate_buckling_analytic.py does), extend that test with 0/90 laminates; tolerance and convergence criteria predeclared
in the docstring. Then T8: compute_prestress_field vs a closed-form centrifugal stress, include a case with non-zero B kappa (offset
reference surface or unsymmetric laminate). One anchor per session if the context grows. lambda1 vs CCX is closed as a diagnostic.


## Handoff 2026-10-19 (12) - START HERE: hypothesis (iv) TE_FLAT slivers, with the localised target of L35
Supersedes (11). Read L34-L35 only. Step A: full suite SLURM 11619070 (sequana_cpu, e4473ee) was still RUNNING at close; read
logs/check/full-slurm-11619070.log tail vs scripts/known-red-full.txt (untracked); outside-registry failures need a 54a6386 worktree
rerun of that test alone; no --record unless justified. Step B is closed (procedure validated, excess localised). Next (iv): split OUR
energy on the TE_FLAT slivers into membrane / drill / transverse shear / bending (element-level, e.g. ours-elem with the drill term
switched off or K_e of element 1025 decomposed), and test the sliver directly: lambda1 with the TE_FLAT strip merged/removed near the
tip (both codes) - the L29 bracket 2.85-3.22 is the prior. Tools: logs/check/tipmodes/t7-v/kselftest.py (stages gen/ccx/ours/split/
elem/local), klocal-dev.sbatch. sequana_cpu_dev accepts only ONE submitted job (maxsubmit 1): no arrays there.

## Handoff 2026-10-19 (11) - START HERE: read full suite 11619070, then localize the K-energy excess
Supersedes (10). Read L32 only. Step A: read tail of logs/check/full-slurm-11619070.log (check.sh now prints ERROR lines and fails on them); compare with scripts/known-red-full.txt (untracked); any failure outside it needs a 54a6386 worktree rerun of that test alone; no --record unless the run justifies it (the stale entry test_blade_modal_gap_converges_with_mesh is one reason).
Step B next (see L33: the 3-8x is a suspect measurement, validate the procedure FIRST with our own eigenmode and a smooth field in both codes): then localize any excess. Add `*EL PRINT, ELSET=T6ALL / ENER` to the kenergy deck (per-element strain energy in CCX) and compute our per-element energy (Rust element K x condensed phi), then compare by element set (TE_FLAT, TE_REINF, TE_PANEL, SW, rest) and by z/L. Also check whether CCX phi translations are mutually consistent at midside nodes. Then hypotheses (iii) MITC4 vs S8R coupling/shear path, (iv) TE_FLAT strip, (v) D E2. SLURM rules unchanged; TMPDIR is $SCRATCH/tmp (check.sh sets its own).

## Handoff 2026-10-19 (10) - START HERE: read the full-suite result, then lambda1 vs CCX
Supersedes (9). Step 1 of (9) is DONE (L31). State: fix 41e9a16 + docs, branch feat/rotating-local-tip-modes, no push; full suite SLURM 11619021 (sequana_cpu, --ntasks=4 only, log logs/check/full-slurm-11619021.log, snapshot logs/check/git-snapshot-20261019T153609, revision 9fb7385 clean; 11619005 and 11619020 were cancelled: 48 cpus requested / bad snapshot export).
Read only L27-L31 and this handoff; mem_search "k-g-prestress-coupling". Rebuild the .so first only if crates/ changed.
Step A - read the suite: tail of the log, then compare against scripts/known-red-full.txt (untracked, not from this feature; do NOT --record unless a real run justifies it, and say which run). Report every failure outside the four documented-red tests; for each, was it red before 41e9a16? (baseline: docs/origin_main_integration_2026-09-30.md and the 54a6386 behaviour; if needed re-run only that test on a worktree of 54a6386 - never infer). Expect K_G-dependent anchors with unsymmetric laminates to move: report numbers, never widen a band.
Step B - the user's priority: lambda1 3.22 vs CCX 2.29 (+41 %), also 2nd/3rd 4.07/4.30 vs 2.57/4.18. Known: prestress now agrees (L27), K_G recipe consistent, load-work ratio MITC4/CCX = 0.844 (L-T7 static), ours stiffer. Candidate splits to test, in order: (i) T7(iv) parked TE-panel K vs CCX (isolates K from K_G); (ii) lambda with K_G held fixed and K swapped / or CCX K_G and our K, on the same mesh; (iii) MITC4 vs S8R membrane-bending coupling (B) and shear-correction path in the sandwich TE panels (T7(iii) xfails: A55 = sum G_k t_k); (iv) TE_FLAT sliver strip sensitivity (L28: 2.85 vs 3.22 bracket); (v) D E2 term deliberately excluded (L28): quantify lambda with it included only on non-sliver elements. Small SLURM jobs, 4 cpus, dev for <20 min (rule below).
Step C (after A/B): T5 closure note, S-4 docstring only if a number moves, TE_FLAT decision, T7(ii), old-mode Rayleigh flag 9.37 vs 7.29 (L29). S-4 criterion: user decided 5 % is his yardstick and 5-6 % between different codes is acceptable context; ElastoDyn is the only arbiter that rotates (BeamDyn 7.56 rpm failed 4/4); do not change the band without him.
SLURM rules (user, 2026-10-19; in AGENTS.md): request ONLY --ntasks=N (no nodes, cpus-per-task or mem); <20 min -> sequana_cpu_dev; longer -> the partition (sequana_cpu or bigmem) with fewer of my running jobs (cap ~4 per partition); several small jobs. Account has a submit/slot cap (AssocMaxSubmitJobLimit): chain with scripts like logs/check/kgrevalid/submit_rest.sh.

## Handoff 2026-10-19 (9) - superseded by (10): targeted K_G re-validation, then full suite
Supersedes (8). State: fix 41e9a16 (K_G prestress (A eps + B kappa)/h, D E2 excluded), doc 32399bb, branch
feat/rotating-local-tip-modes, no push. Validated on the blade production path (L29, SLURM 11618946:
lambda1 3.2244 vs diagnostic 3.2244, B=0 fields 1e-16). Read L27-L30 only; mem_search "k-g-prestress-coupling".
Rebuild first if the .so may be stale: cd crates/aeroelast-py && maturin develop --release (aeroenv).

Step 1 - the tests that led to this diagnosis (K_G-dependent), run as ONE SLURM job, log to logs/check/:
- direct: tests/test_kg_prestress_coupling.py, tests/validation/element/test_plate_buckling_analytic.py,
  tests/validation/element/test_stress_stiffened_solver.py, tests/test_t7_panel_stiffness_parity.py,
  tests/test_iea15mw_s4_rotating_modal.py (the anchor that opened this feature: S5 1st edge red 5.9-6.1 %
  vs 5 %; recorded S-4 at es 0.25: 1st flap -1.4 %, 1st edge -5.6 %, 2nd flap -1.7 %),
  tests/test_iea15mw_s0_omega_zero_consistency.py.
- rotating K_G users: tests/test_iea15mw_sx_campbell.py, tests/test_iea15mw_sx_campbell_cross.py,
  tests/test_rotor_centrifugal_branch.py, tests/test_omega_rebuild_predicate.py,
  tests/validation/rotor/test_rotor_rust_parity.py, tests/test_inertial_rotor_rust_parity.py,
  tests/test_inertial_corotational_eigen_parity.py, tests/validation/blade/test_blade_deloading_vs_reference.py.
- S-5/S-6 (tests/test_iea15mw_s5_oneway.py, s6) only if time allows in the same job; else step 2.
Report every anchor number that moved vs the recorded values (docs/origin_main_integration_2026-09-30.md,
docs/validation_closures.md, test docstrings); never widen a band; a moved anchor is a finding.
Also re-run the local-mode diagnosis that started this (S1-S4): update tools/diagnose_rotating_prestress.py
(:111-163) to asm._rust.compute_prestress_field, then rerun its solve/buckle sbatch at es 0.25 oml to NEW
output names (do not overwrite *-oml-025.npz baselines) and check whether the low-beam-fraction tip modes at
7.56 rpm (z/L 0.94-0.97) disappear.
Step 2 - full suite via SLURM: scripts/check.sh full, compare with scripts/known-red-full.txt (untracked,
unrelated to this feature; re-derive with --record only from a real run and say which).
Step 3 - T5 closure (record the finding, S-4 docstring if a number moves), T7 vs CCX 2.29 residual (operator
x0.71, W ratio 0.844), TE_FLAT sliver decision (L28), T7(ii). Open flag: old-mode Rayleigh 9.37 vs 7.29 (L29).

## Handoff 2026-10-19 (8) - fix implemented and validated on the blade; full re-validation next
Supersedes (7). Read L27-L29. Next session, in order: (1) full suite via SLURM (scripts/check.sh full),
compare with known-red-full; (2) S-4 rotating modal vs OpenFAST v5 MBC3 (K_G changed for unsymmetric
laminates: report every moved number, do not widen); (3) T7 tip buckling vs CCX 2.29 - remaining factor
~0.71 is the operator difference (W ratio 0.844), separate question; (4) update
tools/diagnose_rotating_prestress.py to compute_prestress_field; (5) TE_FLAT sliver strip modelling decision
(L28); T7(ii), T5 still open.

## Handoff 2026-10-19 (7) - cause demonstrated; production fix awaits user decision
Supersedes (6). Read L26-L27. Physical cause demonstrated (L27). Next, only after explicit user decision:
consistent membrane resultant in the K_G stress recovery (mitc4.rs compute_element_stress path used by
assemble_geometric_k_from_disp / assemble_geometric_stiffness; MITC3 equivalent), maturin rebuild, then the
anchors (S-4 rotating modal vs OpenFAST, T7 buckling vs CCX 2.29, full suite). Before that, optional
diagnostic: locate the z/L 0.79 consistent mode (TE_FLAT D E2). T7(ii), T5 open.

## Handoff 2026-10-19 (6) - B*kappa necessary, not sufficient; next: our consistent resultant
Supersedes (5). Read this, L26, and mem_search full-prestress-tensor. Established in L26 (numbers there).
Next, in order, still diagnostic-only (no formulation change until shown):
1. Our own consistent membrane resultant N = A*eps_m + B*kappa + cb*E2 (mitc4.rs W blocks, :1262) from OUR
   static u (kappa=(2/h)e^b1, E2=(4/h^2)e^b2); K-fixed swap with it. Prediction to test: lambda well above
   0.86 only if our kappa also differs from CCX's or eps_m is coupled; report whatever it gives.
2. Explain the remaining membrane-strain difference in the TE_REINF (ours A*eps vs CCX A*eps0) - compare
   eps0 directly, not resultants; suspects: thin unsymmetric TE_REINF stiffness/coupling in K, S8R vs MITC4.
T7(ii) and T5 remain open. Artifacts: t7-iv/{bkappa_ccx.py,bkappa-ccx.npz,bkappa-swap-dev.sbatch,
bkappa-swap-11618888/result.json, bkappa-swap-11618888.log} (ignored, preserve).

## Handoff 2026-10-19 (5) - STOPPED by user; analyze stress-swap 11618827 next
Supersedes (4). One objective next session: identify the physical cause of the 0.52 vs 2.29 buckling gap,
not add tests or fix section shear speculatively. Branch feat/rotating-local-tip-modes; diagnostic commit af7b8fe;
feature base local main 29897f1. Delivery single-pr; no publishing authorized.

**Established physical evidence** (L20, L23):
- Parked 09_46 patch under one identical discrete load: W_MITC4/W_CCX = 0.843789 (MITC4 less compliant).
  This does not clear K for the blade buckling shape; no physical parity band exists.
- The FULL prestress is not the same. Transverse mean N ours/CCX in N/m: panel 46 -430.12/-106.86;
  panel 47 -418.63/+12.12. Matching only N_zz hid this difference. Frames, element mapping and stored
  N_zz reconstruction were checked (errors 0 m, 1.21e-16 and 3.76e-8 respectively).
- Recovery feeding our K_G uses A*epsilon/h, omitting B*kappa from the full unsymmetric laminate
  resultant. This is a candidate implementation defect, NOT yet the demonstrated cause of the gap.
- Sandwich section-shear discrepancy is separate; it does not establish the cause in these no-core tip
  panels. Do not modify production shear or any acceptance band.

**First action**: read this handoff plus L22-L25, recall `full-prestress-tensor`, then inspect
`logs/check/tipmodes/t7-iv/stress-swap-11618827/result.json` and its stdout. The job exited 0 but its values
were NOT read before closing. Same K and free DOFs in all cases:
1. `ours`: confirm baseline lambda ~0.5176621298 and acceptable residuals.
2. `ccx_all`: quantify change with the complete CCX tangent prestress.
3. `ccx_transverse_only`: quantify the effect of replacing just the transverse component globally.
4. `ccx_transverse_tip_only`: same replacement only in 09_42..49_LP_TE_PANEL.
Compare the original-mode destabilizing-energy ratios and Rayleigh factors as well as the new lowest
positive factors; mode switching is not proof of a mechanism. If stress substitution explains the gap,
next isolate why the stress differs (quantify B*kappa, rather than assuming it). If not, examine remaining
stress/geometric-stiffness effects. T7(ii) orthotropic anchor and T5 closure are still open.

**Numerical history / limits**: 11618824 failed at positive-factor selection after reproducing the baseline;
the corrected job uses largest-real inverse eigenvalues, not largest-magnitude. No matrices/stresses were
changed by that correction. Panel job 11618805 failed JSON serialization, fixed with a realistic regression;
11618812 produced the static result. Diagnostic checks: 7 passed, Ruff clean, quick 75 passed/1 xfailed;
existing sandwich 2 strict xfails remain. No Rust rebuild/full suite (production unchanged). Native assessment
was unavailable due undeclared untracked paths; no review approval claimed. scripts/known-red-full.txt is
unrelated and untouched. Scientific logs/output remain local ignored artifacts; preserve them.
Engram stores this active handoff and discoveries; full historical-document mirror remains pending.

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
