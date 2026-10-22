# K16 G1 output-contract preparation
Prepare a non-overwriting full-field input candidate; do not execute BeamDyn or claim G1–G4 closure.
Branch: `audit/b3-s2-beamdyn-static`; baseline: `5f8222f`; delivery: ask-on-risk, local commits only; forecast ~350 authored lines; checks: bootstrapped focused pytest and quick smoke.

## Specs
- S1: "Objetivo: cerrar la viabilidad G1-G4 antes de ejecutar el discriminante." Current bounded unit prepares the G1 candidate; authoritative criterion: `odd/tasks/validation-l3-b3-audit.md`, K16 section only, starting at line 1840. No numerical attribution follows from input preparation.
- S2: "El deck inspeccionado no declara TDzr en todos los nodos. Faltan equivalencia geométrica/de cargas y cotas independientes de incertidumbre. Residual no equivale a GA. No se leyeron nuevos resultados físicos." Check the candidate requests all six nodal TD/RD components; missing reference geometry/directors, provenance equivalence and G2–G4 proof remain explicitly missing, never synthesized or passed.
- S3: "K14 VOID; K15 condicional; LC3 undecided; D1-D3 pendientes." "No tocar contract.yaml ni ejecutar E/F. No repetir K13 completo." Checks: no solves, production changes or verdict changes.
- S4: "Yaml oficial; malla 0.15–0.25 m, 0.125 solo convergencia." "Python con scripts/aeroenv.sh y -u." "SLURM solo --ntasks=N; M1 0.125 en sequana_cpu, 16 tasks, nunca dev." "Esperas background sin polling. Sin push." No mesh or SLURM execution is needed for this text-only unit.
- S5: "Preservar todos los directorios untracked de referencias BeamDyn." Check explicit isolated destination, refusal to overwrite, unchanged input references and no binary execution. Existing build/run/post behavior must remain unchanged.
- S6: "bueno dale arranca". Authorization to start prerequisite implementation, not permission to run the physical discriminant. Assumption: an additive, generate-only candidate tool is the smallest safe implementation; no modification to the existing generator is required.

## Tasks
- T1 | S1–S6 | delegated writer (context-backstop trigger) | done | Additive generate-only tool/help + 11 public CLI tests; observed RED/GREEN, Ruff and smoke pass; local work-unit commit pending.
- T2 | S1–S6 | parent reconciliation | in_progress | Reconcile exact diff, final spot-check and G1–G4 missing proof; record local work-unit commit; no physical closure claim.

## Log
### L1 — Original request (verbatim)
```text
Retomamos B3-3 en audit/b3-s2-beamdyn-static, HEAD 5f8222f, local sin push.

   Leé mem_context y hacé una búsqueda acotada "K16" (observaciones 68 y 69;
   handoff 70). En odd/tasks/validation-l3-b3-audit.md localizá con grep y
   leé SOLO la sección K16 del final; recurrí a K15 RESULT y al Handoff
   2026-10-22 únicamente si hace falta. No releas el documento completo.

   K16 preregistrado en 586bd16; smoke documentado en 5f8222f:
   75 passed, 1 xfailed, exit 0.
   Log logs/check/quick-20261022T093443.log. No es cierre físico.

   Objetivo: cerrar la viabilidad G1-G4 antes de ejecutar el discriminante.
   El deck inspeccionado no declara TDzr en todos los nodos. Faltan
   equivalencia geométrica/de cargas y cotas independientes de incertidumbre.
   Residual no equivale a GA. No se leyeron nuevos resultados físicos.

   IMPORTANTE: las referencias a openfast-src fuera del repo son provisionales,
   no evidencia admitida: el subagente lo consultó sin autorización.
   Pedí autorización explícita antes de volver a inspeccionar esa ruta.

   K14 VOID; K15 condicional; LC3 undecided; D1-D3 pendientes.
   No tocar contract.yaml ni ejecutar E/F. No repetir K13 completo.
   Yaml oficial; malla 0.15–0.25 m, 0.125 solo convergencia.
   Python con scripts/aeroenv.sh y -u.
   SLURM solo --ntasks=N; M1 0.125 en sequana_cpu, 16 tasks, nunca dev.
   Esperas background sin polling. Sin push.
   Preservar todos los directorios untracked de referencias BeamDyn.
```
### L2 — Subsequent user instructions (verbatim, chronological)
```text
si claro,
si
instala openfast con conda, estoy seguro ue ya eso esta instalado en el env que estamos usando, o que se puede instalar
bueno dale arranca
```
The first two replies granted the explicitly presented read-only installation and recipe paths. Source-root read access was separately granted through the choice UI. No outside write/execute authority was granted. See Engram 71–77 for exact target bindings and provenance evidence.

### L3 — Reconciled baseline and implementation boundary
- Source-only review supports all-node TD/RD, WM rotations and Gaussian FE/GLL motion nodes; no configuration-only export of undeformed positions/directors was established. The candidate must state these gaps, not invent channels.
- Installed BeamDyn matches local conda-forge `openfast-5.0.0-h9ad4a36_1` manifests. Recipe declares upstream v5.0.0 tarball; tarball-to-checkout byte equivalence remains unverified. No reinstall required. Main AeroElast Python uses a separate venv.
- Existing helpers live in `tools/bd_static_reference.py`; exploration mapped generation at lines 106–124, exclusive destination creation, and text helpers. The new tool must reuse suitable pure helpers without changing existing commands or their shared code.
- Route: bounded writer because the runtime signaled the parent context backstop. Allowed edits are the new tool and its new tests only; the parent owns this document, reconciliation and commits.
- Functional checks are text generation/CLI tests, not physics. The quick smoke (~13 s nominal; prior run ~103 s wall time) is sufficient repository smoke for this unit. No full physical suite or binary invocation.

### L4 — Observed checks and corrected reference-export feasibility
- Writer T1 RED: 10 failed, exit 1 (new CLI absent); GREEN: 10 passed in 4.88 s, exit 0. One assertion-whitespace correction; Ruff initially rejected `RUFF_NO_CACHE=1`, corrected to `true`, final exit 0. Parent spot-check: 10 passed in 4.52 s, exit 0.
- Quick smoke job-1: exit 0 after 17 s; 75 passed, 1 xfailed, 699 warnings in 13.93 s; Ruff clean and mesh-size OK. Wrapper `logs/check/k16-g1-candidate-quick.log`, canonical `logs/check/quick-20261022T101326.log`. Smoke is not physical validation.
- Native assess over untracked additions was unassessable (explicit untracked scope required); returned plan keeps writer self-verification, no independent verifier without a named high-risk item. Reassess the exact staged work unit, not the accumulated branch.
- Authorized source-only follow-up corrects the earlier incomplete export assessment: `BeamDyn_IO.f90:2049–2060` writes `Init_Nodes_E<n>` with initial XYZ and WM orientations when SumPrint is enabled. Parent spot-check confirms per-element arrays and missing global indexing. `ES18.5` formatting is finite-precision and needs an independently propagated serialization bound before physical use.
- The existing candidate already sets SumPrint=True and tests all six TD/RD components, including RDxr. Add source-supported summary schema metadata, distinguishing available-by-design from acquired/validated. Do not implement a physical-output parser or consume outputs in this unit. No source instrumentation appears necessary to request the initial geometry; ordering, units, coverage and conversions remain unverified.
- Metadata follow-up RED: 1 failed, 10 passed in 4.92 s, exit 1; GREEN: 11 passed in 6.18 s, exit 0; Ruff exit 0. Parent final spot-check: 11 passed in 7.14 s, exit 0. Existing commands/helpers unchanged; differential does not apply. The earlier smoke precedes this metadata/help-only addition; no repeated smoke.
- Public CLI generated two persistent, non-overwriting input-only candidates successfully (exit 0): `logs/check/k16-g1-inputs-LC3-p60-s0.0001` (edge) and `logs/check/k16-g1-inputs-LC2-p60-s0.0001` (flap control). Each has beam.dat, blade.dat, driver.dvr, candidate.json; neither has been executed. These are not physical evidence or matched-mechanics references.
