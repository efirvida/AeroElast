# CCX FSI parity — parked V50 (issue #10)

**Goal**: run the failing parked V50 FSI case twice against the *same* BEM fluid
participant — once with the aeroelast solid, once with CalculiX+preCICE
(`ccx_preCICE`) — and compare tip-displacement histories. Verdict rule: if CCX
settles (~4.7–8.3 m reference, `docs/validation_data/generated/parked_flag_v50_response_summary.csv`)
where aeroelast diverges (2168 m, issue #10), the defect is in our solid/coupling;
if both diverge alike, it is the aero/coupling configuration.

**Context / evidence**
- Case files: `tests/IEA15MW/parked_v50/{solid_v50.yaml,fluid_v50.yaml,precice-config.xml}`
  (StressStiffenedDynamicFSI, V50/pitch90/yaw8/omega0, dt_window 1e-2, max-iter 70,
  IQN-ILS, force_ramp 2.0 s, Rayleigh zeta 0.03).
- Failing run: `$SCRATCH/parked_v50_rerun/case` (preCICE 1519 windows, it/win 9.4,
  tip |Y| 2168 m at t≈10 s).
- Fluid is the BEM preCICE participant (`fem-shell-bem-fsi`) over sockets — it is
  reused unchanged; only the solid participant is swapped.
- `ccx_preCICE` already installed in venv (Makefile.sdumont `calculix` target);
  runtime needs `LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64:$VENV/lib:$VENV/lib64:/usr/lib64/psm2-compat`
  (Engram `env/ccx-precice-adapter-venv`).
- CCX writer: `aeroelast.core.mesh.io.writers.write_ccx_mesh` (composite layup,
  node sets, `*DYNAMIC`, dt/t_end, shell element type); helpers `tests/_ccx_io.py`.
- Adapter (`precice/calculix-adapter`, cloned at
  `/tmp/pi-github-repos/runtime-VfFMZy/82686281b7024400b152d85a31de3ea5c04c331911b2e0b4d79aa1eee188ce5f`):
  `config.yml` with `nodes-mesh`, patch → `*NSET` (N-prefix added automatically),
  read-data `Forces` (`*CLOAD`), write-data `Displacements`/`Velocities`;
  data names are fixed keywords and must match the XML; all interface nodes must
  carry a zero `*CLOAD` at step start; `*DYNAMIC` is mandatory for FSI.

**Design decisions (apples-to-apples)**
- Two coupled runs share the same BEM fluid config; only the solid participant
  (and its data names in the XML) changes.
- CCX data keywords are plural (`Forces`/`Displacements`); our XML uses singular
  (`Force`/`Displacement`) → parity XML variant; BEM side must accept the renamed
  data (config or minimal patch — task 1 resolves where the names live).
- `Velocity` exchange: solid_v50 writes it, fluid_v50 has `velocity_data: null`;
  XML still carries the data — reconcile what actually ran, then make both parity
  runs consistent (prefer: no Velocity if BEM tolerates it).
- Force ramp 2.0 s is applied solid-side in aeroelast; CCX cannot ramp → apply
  ramp fluid-side if trivially possible, else disable in both parity runs
  (divergence grows after t=2 s, ramp only shapes startup).
- Rayleigh 3 %: verify CCX 2.20 damping support (`*DAMPING`/HHT α); fallback =
  disable damping in both parity runs, physics stays identical across solids.
- Run window: ~15–20 s is enough (aeroelast diverges by t≈10 s; reference
  settles by t≈5 s).

**Task 1 — resolved (reconciliation)**
- BEM preCICE data names are YAML-configurable
  (`solvers/bem/fsi_participant.py:1473-1482`: `force_data`, `displacement_data`,
  `velocity_data`; `velocity_data: null` already used by `fluid_v50.yaml`, guard at
  line 492) → **no code patch**: parity fluid variant sets
  `force_data: "Forces"`, `displacement_data: "Displacements"`.
- XML drift confirmed: `precice-config.xml` (Sep 9) still declares/exchanges
  `Velocity` while `fluid_v50.yaml` disables velocity coupling and
  `solid_v50.yaml` writes only `Displacement`; the failing run dir kept no XML
  copy (`parked_v50_rerun/{case,precice-run}`). → parity XML variant: data renamed
  to `Forces`/`Displacements`, Velocity data/exchange removed.
- Rayleigh damping: `*DAMPING` (α/β) is compiled into `ccx_preCICE` (rayleigh_*
  symbols present; "no Rayleigh mass damping allowed" belongs to another
  procedure) → replicate aeroelast's α/β for zeta 0.03 instead of disabling.
- Force ramp 2.0 s: solid-side only in aeroelast → **off in both parity runs**
  (`force_ramp_time: 0` in the aeroelast variant), documented; the historical
  failing run kept 2.0 s and still diverged after t=2 s.
- Verdict setup: parity campaign = run A (aeroelast, parity XML+config) and
  run B (CCX, same XML, same fluid); the historical failing run remains context.

**Tasks**
- [x] 1. Reconcile the case as-run: where BEM data names live, actual XML in
      `parked_v50_rerun`, Velocity/ramp/damping decisions recorded here.
      → resolved above; campaign = two coupled runs with a parity XML variant.
- [x] 2. CCX deck exporter → `tools/ccx_fsi_parity_parked_v50.py` (deck,
      `config.yml`, parity XML, both parity YAMLs) + the `DynamicFSI` branch in
      `write_ccx_mesh` (`tests/test_ccx_writer_fsi_deck.py`). Verified: the deck
      runs in real CCX (`Job finished`, FRD written) and the Rayleigh mapping
      reproduces zeta 0.03 in both damping modes.
- [x] 3. Cheap structural gate → `tests/test_ccx_fsi_parity_parked_v50.py`
      (5 passed, ~106 s): artifacts parse, adapter cards present, damping
      mapping, modal CCX-vs-AeroElast on the same mesh, and a real CCX run of
      the generated dynamic deck.
- [x] 4. Coupled-run wiring → `tests/IEA15MW/parked_v50/ccx_parity/run_cost_measure.srm`
      (deadline-aware SLURM script: generates the 0.25 m deck, times the
      CCX-only increments, then runs the real coupled pair with the solid
      acceptor first) plus the generator options it needs: `--alpha/--beta`
      (skip the 290k-DOF modal solve), `--abs-paths` (run dir outside the repo)
      and a self-contained `exchange-directory="."`.
- [ ] 5. Short coupled run on SLURM + preCICE monitoring. Cost measurement
      submitted as job 11606103 on `sequana_cpu_dev` (20 min): production deck
      at 0.25 m generates in 156 s → 32 336 nodes / 33 473 S8R elements
      (~194k DOF). Measures CCX-only seconds per increment and then the real
      coupled seconds per window to project the full campaign.
- [ ] 6. Comparison + verdict: tip-disp extraction for both runs, metrics
      (peak / settled mean / diverges), pytest that skips when outputs are
      absent, results doc, evidence appended to issue #10.

**Findings and incidents (task 2/3)**
- Delegate to the official preCICE calculus-adapter requires the interface nodes
  loaded in every spatial direction at step start; the writer now emits the
  three zero `*CLOAD` lines per interface set (`NALLOUTERSHELLNODS, 1/2/3, 0.0`).
- CCX accepts `*DAMPING, ALPHA=<mass>, BETA=<stiffness>`; verified numerically
  that the mapped coefficients reproduce the case's 3 % in modes 1 and 2. The
  aeroelast solver names them the other way around internally (its `alpha` is
  the stiffness term), so the mapping is explicit in the tool.
- Element choice matters: `quadratic=False` (S4R + smeared equivalent laminate)
  gave 150 % modal gaps against MITC4Composite at the gate mesh; S8R with a real
  `*SHELL SECTION, COMPOSITE` (the configuration the repo's blade parity tests
  already validate) fixed it. The deck uses S8R.
- `_aeroelast.compute_rayleigh_auto` (Rust fast path) fails on this host with
  `PETSc error 95` in `EPSSolve`/`EPSGetEigenpair` although `modal_solve_coo`
  solves the same problem on the same matrices. The tool falls back to the
  solver's closed-form formulas. Candidate upstream defect, not blocking.
- Environment: `python` on PATH resolves to the pi-lens pip-tools interpreter,
  not the project venv; the venv also needs the GCC 14 libstdc++ in
  `LD_LIBRARY_PATH` (`CXXABI_1.3.15`). Run parity commands with
  `$VENV/bin/python` and that `LD_LIBRARY_PATH`.

**Acceptance**: both coupled runs complete the window; comparison report shows
either "CCX settles where aeroelast diverges" (solid-side bug, narrows issue #10)
or "both diverge" (aero/config problem); all artifacts reproducible from the repo.

**Risks**
- Degenerate blade element (CCX `nonpositive jacobian`, element 2790 class) —
  mesh at 0.25 m passed CCX modal before; if the deck hits it, regenerate with
  the workaround from `docs/origin_main_integration_2026-09-30.md`.
- CCX implicit dynamics is HHT (α-method), not Newmark-β — small transient
  differences acceptable; verdict rests on settle-vs-diverge, not pointwise match.
- MITC4Composite vs S4R laminate equivalence already validated by existing
  static/modal CCX parity tests (0.4–4 % gaps).

## Cost measurement — job 11606103 (`sequana_cpu_dev`, 20 min)

Xeon Gold 6252, 8 CPUs, `OMP_NUM_THREADS=8`.

| Quantity | Measured |
| --- | --- |
| Production deck generation (0.25 m, Rayleigh override) | 156 s |
| Deck mesh | 32 336 nodes / 33 473 S8R elements |
| Coupling initialisation | **works**: preCICE logged `TimeWindow 1, Iteration 1` and the Flap-Tip watch-point |
| Windows completed at 0.25 m in 521 s | **0** |
| CCX-only progress in 420 s at 0.25 m | first increment still running |
| CCX-only peak RSS | **11.4 GB** |
| Local scaling 1.0 m (3 043 nodes) | 13 increments in 600 s (~46 s/increment) |
| Local scaling 0.5 m (9 277 nodes) | 7 increments in 601 s (~86 s/increment) |
| CPU inside CCX | ~152 %; "Using 1 cpu for spooles" |

Two problems, cost first:

1. **Cost**: one coupling window is ~9 preCICE iterations, i.e. ~9 CCX increments.
   At 0.5 m that is already ~13 min per window (~9 h for a 40-window run), and at
   0.25 m the window never closed. A direct solver cannot carry this.
2. **The deck diverges under load**: with a representative 100 N per interface node
   the local runs blow up (largest displacement increment 67 m at 1.0 m, 2.4e4 m at
   0.5 m), and the 0.25 m log reports a residual force of 5.3e13. The modal gate at
   2.0 m matched AeroElast, so the *stiffness* is plausible — the suspicion is the
   composite shell model itself: CalculiX expands `*SHELL SECTION, COMPOSITE` into
   internal layers (node ids reached 647 462 in a 32 336-node mesh, i.e. ~10–20x
   expansion), which both inflates the model and can wreck conditioning.

**Decision taken**: abandon the layered-composite `*SHELL SECTION` deck as the
parity solid. Next candidate is `*SHELL GENERAL SECTION` with the explicit CLT
A/B/D matrices the solver already computes (`Laminate`), which keeps the layup
mechanics in a single layer, removes the expansion, and stays cheap. Verify with
the existing gate (modal vs AeroElast plus a static load at 1.0 m) before spending
another coupled run.

New risks recorded: (a) the coarse-mesh route needs the *fluid* mesh regenerated at
the same element size, so the generator must patch `element_size` in the parity
YAMLs; (b) the divergence must be re-checked on the replacement element, since it
may be a real modelling defect and not just cost.
