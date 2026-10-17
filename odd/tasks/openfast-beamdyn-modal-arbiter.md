# Feature: self-consistent dynamic arbiter for the IEA-15 v1.1 blade (OpenFAST v5 + BeamDyn), runs only
Branch: main. Parent: `odd/tasks/validation-l3-b3-audit.md` ("B3-2 RESULT"). Analysis is a SEPARATE session (maintainer, 2026-10-17).

## Specs
- S1: "si claro tienes mi autorizacion para hacer todas las corridas que entiendas necesarias"
- S2: "pero el analisis lo hacemos en otra sesion de trabajo" (this unit produces runs and raw numbers only; no store edits, no verdicts, no conclusions about the shell).
- S3: from the previous report, accepted by S1: "OpenFAST v5 con BeamDyn sobre el deck v1.1" (blade structure from the BeamDyn 6x6 deck instead of ElastoDyn's mode-shape representation), "S-4 1F/2F y reabrir los veredictos de V-02".

## Tasks
- [ ] T1 locate or obtain an OpenFAST v5.0.0 binary (and a BeamDyn-enabled IEA-15 v1.1 deck); record paths.
- [ ] T2 blade-only modal runs, BeamDyn blade: parked (0 rpm) and 7.56 rpm, collective (MBC3) frequencies for 1F, 1E, 2F, 2E (+ torsion if present).
- [ ] T3 same two cases with the ElastoDyn blade of the same deck as a control (must reproduce the stored S-4 targets 0.544/0.739/1.610 parked, 0.5666/0.7461/1.6372 rated within ~1 %; if not, report before anything else).
- [ ] T4 raw result table + handoff block in this file (paths, job ids, commands), no analysis.

## Log
- L1: "continuamos con la validacion ya terminaros las simulaciones de S4" (session origin); L2: "si claro tienes mi autorizacion para hacer todas las corridas que entiendas necesarias"; L3: "pero el analisis lo hacemos en otra sesion de trabajo"
