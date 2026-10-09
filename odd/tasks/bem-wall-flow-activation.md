# Feature: activate the BEM wall-flow moment realisation in production (#16, item P2 of #18)

Status: T1 in progress (2026-10-09). Nothing committed yet.
Owner: this session (2026-10-09)
Related: issue **#16** (roadmap item `P2` of **#18**), `docs/validation/gaps.yaml` id
`moment_realization_over_delivers`, `docs/validation/rows/31-tube_moment_realization.yaml`
(group 31: the single-cell tube at 0.3091% against Bredt), `tests/test_multicell_shear_flow.py`
(the multi-cell system, pinned against a hand-assembled Bredt-Batho), issue #11 (the
realisation defect, closed) and issue #17 (per-process node ordering, adjacent).
The decision was taken by the user on 2026-10-09: **activate**, then measure.

## Measured on the production path (2026-10-09, `element_size: 0.25`, the campaign's own setting)

Probe `$SCRATCH/bfs16/probe_ring_counts.py` (outside the repo): it calls
`aeroelast.cli.run_bem_fsi._build_mesh` with the campaign's own `mesh` section and then
builds the production `ForceProjector` with the deck's property map.

| quantity | before (`a57f4c5`-era filter) | after T3 (`subset_to_nodes`) |
| --- | --- | --- |
| coupling mesh | 27609 nodes, **0 elements** | 27609 nodes, **27598 elements** |
| full mesh kept as `viz_mesh` | 32325 nodes, 33462 elements | unchanged |
| ring sections usable | **0 of 671** | **671 of 671** |
| `from_element_properties` | n/a (nothing built) | 671 of 671 |
| rings whose walls all resolved `S = G*t` | n/a | 671 of 671 (`0` rings with the uniform `S = 1.0`) |
| cells per ring | n/a | **min 1, max 1** (the coupling set excludes the webs) |
| walls per ring | n/a | min 16, max 68 |
| property coverage of the coupling elements | n/a | **27598 of 27598 (100.00 %)** |
| projector init (once per process) | 0.4 s (no rings) | 50.1 s |

Two numbers matter beyond the closure: **the coverage is total**, so the latent gate defect is
not stepped on (the map resolves every coupling element through its own section set), and the
rings are **single-cell everywhere**, which is the limit the store entry has to declare. The
50.1 s init is paid once per participant process (`ForceProjector.__init__` precomputes
`_strip_ring_sections`; `project()` reuses them), not per window.

## Why this exists

Issue #16 measured that the multi-cell Bredt-Batho realisation **never runs in production**:
on the fluid coupling mesh the projector gets 27609 nodes and **0 elements**, so every
`ring_section` call raises `ValueError: no element chordwise edge joins two ring nodes` and
**0 of 671 ring sections are usable**. Every strip falls back to the minimum-norm
`_distribute` field, which is exactly the field `moment_realization_over_delivers` says
**over-delivers rigid torsion** where a closed section carries wall shear flow.

The closes-when of the item is a decision plus a register that matches production. The user
chose activation. What activation can reach **in production** is the *single-cell skin* wall
flow, not the multi-cell split: the coupling node set (`allOuterShellNods`, used by every
campaign and smoke YAML) excludes the shear-web nodes on purpose, so the ring's closed loop is
the outer skin. That limit is part of the deliverable, not a failure of it.

## The measured state that authorises the change

| Fact | Evidence |
| --- | --- |
| The acople is node-only: `MeshModel(nodes=filtered_nodes)`, elements dropped | `src/aeroelast/cli/run_bem_fsi.py:276-294` |
| The participant builds the projector with no property map, and its own comment says the caller must hand it over for a stiffness-resolved split | `src/aeroelast/solvers/bem/fsi_participant.py:356-370` |
| The gate refuses a ring without properties, so the strip falls back | `src/aeroelast/solvers/bem/force_projection.py:1036-1048`, `:281`, `:328` |
| The generated mesh already carries the skin elements and their section membership | `src/aeroelast/core/mesh/generators.py:1291-1305`; `models/blade/numad/mesh_gen/mesh_gen.py:560-612` (`allOuterShellEls` is the union of the skin section sets, so skin elements keep their own `elementSet` key) |
| Both paths reproduce the same `F` and `M`: the change is the **distribution**, not the resultant | `docs/validation/gaps.yaml` id `moment_realization_over_delivers`: "the realised moment is exact and the net force is zero" |

Because the resultants are identical, the only measurable consequence of activation is the
**spatial distribution** of the applied load, i.e. the elastic twist response. A prediction is
written below so the measurement can falsify it.

### Prediction (falsifiable, the point of T4)

The minimum-norm field over-delivers rigid torsion, so removing it should **reduce** the
over-delivered elastic twist. In the one-way de-loading table the blade currently unloads
**twice** the published reference (-26.36% thrust / -15.74% power against Zhou et al. 2025
Table 6's -13.04% / -8.38%). If the over-delivered twist is a contributor to that gap, the
wall-flow realisation must move those two numbers **toward** -13.04% / -8.38%. If it does not
move them, the load-distribution hypothesis for that gap is dead and the gap stays attributed to
`radii_datum_definition_bias`. Either outcome is a result.

## Design

Four changes, in this order:

1. **`MeshModel` gains a node-subset view that keeps elements.** A new mesh whose nodes are
   exactly the requested ones (in the order the caller gives, so the coupling node list is
   unchanged) plus every *source* element all of whose nodes are in that set, with `NodeSet`
   and `ElementSet` membership restricted the way `extract_submesh` already does. A mesh with no
   surviving element stays nodes-only (the current behaviour).
2. **The CLI filters through it.** `_build_mesh` (`run_bem_fsi.py`) keeps using
   `coupling_node_set` for the node list but takes the elements along.
3. **The properties reach the projector.** The CLI builds the Rust properties map from the very
   generator that produced the mesh (`build_rust_properties(generator.numad_mesh_data)`, with the
   `_blade_N` suffix expansion for `RotorMesh`), and hands it to `build_from_config`, which hands
   it to `BEMFSIParticipant`, which hands it to `ForceProjector`. Default `None` keeps every
   existing direct caller bit-identical.
4. **The measurement.** The coupled gate case on both variants, plus the one-way de-loading
   table as the independent arbiter.

Explicitly *not* in this feature: adding the shear-web nodes to the coupling mesh (that is the
multi-cell path, a much larger physical change nobody has asked for), and fixing the latent gate
defect found below (it gets its own issue).

### Latent defect found while exploring (own issue, not absorbed)

`from_element_properties` is set from `element_properties is not None`
(`force_projection.py:387`), **not** from whether any wall actually resolved a property: a map
whose keys do not match the projection mesh's element sets leaves every wall at `S = 1.0` and the
ring still claims a physical split. `_element_property_by_id` silently `continue`s on a missing
set (`:1209-1222`). Consequence for this feature: the map we pass *must* match the sets the
filtered mesh carries, or the "activated" path is geometric-only while reporting otherwise. That
is a guard condition of T3 and a separate issue at close.

## Tasks

**T1 - RED: the production mesh cannot realise a ring.** Entry: this note; the design above.
Write a test that drives the production builder
(`aeroelast.cli.run_bem_fsi._build_mesh`) on a coarse `BladeMesh` config of the official deck
with `coupling_node_set: allOuterShellNods`, then builds a `ForceProjector` with the deck's
property map and asserts at least one strip's ring section is realisable. Closes when the test
exists and **fails** with the current code, recording in the failure message how many of how many
ring sections were usable (expected 0). Evidence: the failing run, the count.

**T2 - the node-subset view.** `src/aeroelast/core/mesh/model.py` (new method next to
`extract_submesh`) plus its unit test on a synthetic mesh (skin elements + web elements + two
node sets). Closes when the method preserves the caller's node order, keeps exactly the fully
contained elements, restricts both set kinds, and the unit test pins each of those.

**T3 - the CLI carries the elements.** `run_bem_fsi.py`: filter through the new view and build
the property map from the generator. Closes when, on the coarse production config, the coupling
mesh carries elements, the projector built by the CLI path has realisable rings, and T1 turns
GREEN. **Done** (`subset_to_nodes` in the filter, `gen_cfg` init for the pre-existing
"possibly unbound" finding): T1's test GREEN at both sizes, the campaign-scale probe above,
plus a second test pinning the filter's node list, ids and element-set membership.

**T4 - `element_properties` reaches the projector.** `fsi_participant.py` (`__init__` kwarg,
the factory, the comment block) with the default keeping today's behaviour, plus
`_element_properties_for` in the CLI (map restricted to the mesh's element sets, withheld
unless the kept sets cover **every** coupling element) and the 3-tuple return of `_build_mesh`.
Closes when: the factory-built participant on the production mesh realises rings; **`sum F` and
`sum M` of the applied per-strip load are identical between the min-norm and the wall-flow
path**; and a projector without the map still falls back.

**Done 2026-10-09** (coarse production mesh, `element_size: 1.0`, 2828 nodes): 5 tests green in
36.8 s. The invariance guard measured `|dF| = 1.203545e-10 N` (relative **1.403e-16**) and
`|dM| = 3.895602e-08 N.m` (relative **5.984e-16**) - machine precision, seven orders inside the
`1e-9` bound. The *distribution* is what moves: relative L2 **92.22 %**, max `|df| = 5063.7 N`,
and **2828 of 2828** nodes changed. The property-less element-bearing projector is bit-equal
(`np.array_equal`) to the pre-#16 node-only projector, so the fallback is preserved exactly.

**T5 - the measurement.** Run the coupled gate case (`tests/smoke_fix/base_fix/`,
`tests/run_step1b_smoke.srm`) twice: minimum-norm (map withheld) and wall-flow. Report tip twist,
root moments and the torsion response, plus the one-way de-loading table against Zhou et al.
2025 to test the prediction above. Record the moved baseline the way #19/#26 did (two moves
declared). Closes when the numbers are in this note with the command and the revision. No
tolerance is widened and no anchor is retuned to pass; a moved number is the finding.

**T6 - store and closure.** Row(s) for the activated realisation, `gaps.yaml`
`moment_realization_over_delivers` narrowed to what production now does (and the multi-cell limit
it still does not do: web nodes are not in the coupling mesh), the section in
`docs/validation_closures.md`, the comment and close of #16, and the new issue for the latent
gate defect. Closes when the register matches production and #16 is closed.

## Traps

- **Environment.** Every command that imports `aeroelast` needs
  `bash -lc 'module load glu gcc/14.2.0_sequana; export LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64:$LD_LIBRARY_PATH; ...'`
  or the import dies with `CXXABI_1.3.15`. `scripts/aeroenv.sh <cmd>` wraps it.
- **The tree is shared.** There are staged/modified files from another session (`CLAUDE.md ->
  AGENTS.md`, `docs/blade_input_divergence_utd_vs_official.md`,
  `tests/test_iea15mw_v05_structural_properties.py`, `odd/tasks/token-efficiency.md`). Never
  `git add -A`; commit with explicit pathspecs.
- **Mesh size.** `element_size: 0.25` on the official deck gives 27609 coupling nodes; the T1-T4
  tests use a coarser size on purpose. The campaign-scale count (the 0/671 number) is measured
  separately, not inside the test.
- **`MeshElement` ids are a global counter.** New elements built for the view get new ids; the
  property lookup is by element-set membership, never by id, so this is safe - but do not add an
  id-keyed shortcut.
- **`NodeSet.nodes` is a dict built from a `set`**, so the coupling node order in `_build_mesh`
  is not stable across processes. The view keeps the caller's order; do not "fix" the ordering
  here (that is #17's territory) beyond keeping the current list.
- **Never run `coherence`** on the validation store (#24): it rewrites the row files and wipes
  hand-written prose. `tools/tests` runs with `-m "not slow"` until #24 is fixed.

## Commit plan (work units)

One commit per task, Conventional Commits, branch `integrate/origin-main-2026-09-30`:
`test(mesh): ...` (T1), `feat(mesh): ...` (T2), `feat(bem): ...` (T3+T4 or split), `docs(odd):`
(T5), `docs(validation)+chore(store):` (T6, with the publication).
