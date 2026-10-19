# rotating-local-tip-modes

Why do low-beam-fraction (local) modes appear at 7.56 rpm but not parked, with the residual at
z/L 0.94-0.97? IEA-15 official blade, es 0.25, oml and mid. Follows `odd/tasks/s4-mode-classifier.md`.

## Specs

- S1. "Abrir los `.npz`, ubicar los nodos con residuo y ver si coinciden con el cambio de tamaño de malla en la punta o con un borde libre."
- S2. "Recuperar las tensiones de membrana centrífugas en esa zona y compararlas con la referencia analítica (σ ≈ ρω²∫r dA / A). Si da compresión o un pico donde debería haber tracción suave, el problema está en K_G."
- S3. "Repetir el modal con la rotación prendida y K_G apagado, y con K_SP apagado, para aislar cuál de los dos genera los modos."
- S4. "Mirar la tendencia al refinar la malla: 0.25 contra 0.15 m. Un artefacto numérico se mueve con la malla; un modo físico converge."
- S5. "El rojo de S-4 en el 1er modo de borde (5.9–6.1 % contra un límite de 5 %) es un hallazgo de precisión. No se toca la banda."

## Tasks

- T1 [S1] locate local-mode residual (elements, section, skin/web, tip refinement stations) - inline - commit: (tools, this handoff)
- T2 [S2] centrifugal prestress at 7.56 rpm: per-element membrane resultants N vs the section tension T(r) - inline + sbatch - commit: -
- T3 [S3] rated modal with K_G off, and with K_SP off - inline + sbatch - commit: -
- T4 [S4] es 0.25 vs 0.15 local-mode frequencies and beam fractions - sbatch - commit: -
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

## Handoff 2026-10-18 - jobs submitted, nothing analysed yet
SLURM runs on after this session closes. Next session:
1. `squeue -u $USER`; for each job grep its log (`logs/check/tipmodes/*-<jobid>.log`) for
   `prestress solve|wrote|Traceback|SystemExit`. A `SystemExit` about mesh/free-DOF mismatch means
   the regenerated mesh differs from the stored modes: re-run `diagnose_mode_beam_projection.py solve`.
2. T1/T2: `scripts/aeroenv.sh python tools/diagnose_rotating_prestress.py analyse logs/check/tipmodes/prestress-oml-025.npz`
   (then mid). Sanity first: `rq` must equal `w2` per mode (same K_G as the modal), `T_fe/T_ref` ~ 1
   inboard. Then: is the KG term of the local modes large and negative, and are its top elements in
   the tip stations z/L 0.94-0.97 with Nmin < 0? Which set (skin/web), h, area (slivers)?
3. T3: `diagnose_mode_beam_projection.py analyse logs/check/tipmodes/modes-oml-025-nokg.npz` (and
   noksp, mid nokg). Local modes gone without K_G -> K_G confirmed.
4. T4: 11618441 (es 0.15, up to 4 h) -> analyse `modes-oml-015.npz`; compare local-mode freqs and
   beam fractions with es 0.25 (artefact moves with mesh, physics converges).
Open hypothesis to test if K_G is confirmed: recovered sigma_m = (A/h) eps_m drops the B*kappa
membrane-bending coupling (large for oml-offset laminates), and/or a real compressive N from
prebend bending near the tip.
