# Feature: blade skin on the laminate mid-surface (option) + S-1 EA/mass/edge diagnosis
Branch: main (maintainer works on main; commits not pushed). Parent audit: `odd/tasks/validation-l3-b3-audit.md` ("B3-2 S-1 diagnosis", commit 28c7576).

## Specs
- S1: "continua con todo tiees mi aprovacion para hacer todo" (the approval covers the mid-surface mesh option and the EA/mass/edge diagnoses proposed in the previous turn).
- S2: "hacerlo con una opción de `BladeMesh` apagada por defecto" (default behaviour of `BladeMesh` stays the OML mesh; the OML mesh remains for FSI).
- S3: "re-mediría S-1, S-2, S-4 y V-03 uno por uno, con el mismo comando antes y después, y reportaría lo que se mueva" (no band is widened; a moved anchor is reported as a finding).
- S4: "Un calculador independiente del layup del yaml, sin malla ni código de aeroelast, para ver si el exceso de 19 % de EA ya está en el yaml." 
- S5: "Una comparación del EI edge usando el principal de BeamDyn con el twist ya rotado."
- S6: maintainer rules (previous handoff): every L3 tolerance <= 5 %; failure means find the error or a physical explanation, never widen; element_size 0.5 only as the coarse point of a study; commit always.

## Tasks
- [ ] T1 (S2) `BladeMesh` option `skin_reference` ("oml" default | "mid"), shared inset function, unit tests (default unchanged, mid moves skin inward, no inverted element). Route: worker.
- [ ] T2 (S3) before/after S-1, S-2, S-4, V-03 (es 0.25) with the option off/on; report. Route: parent via sbatch.
- [ ] T3 (S4) `tools/` independent yaml-layup calculator of mass and EA vs 6x6 and vs the mesh CLT sum. Route: worker.
- [ ] T4 (S5) EI edge vs BeamDyn principal-axes-rotated reference. Route: parent/worker.
- [ ] T5 close: store rows/docs, AGENTS.md S-4 row (1E -5.5 %), handoff.

## Log
- L1: "continua con todo tiees mi aprovacion para hacer todo"
