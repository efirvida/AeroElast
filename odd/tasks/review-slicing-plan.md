# Review slicing plan

The branch `test/physical-correctness` is 159 commits and 154 files against `origin/main`
(`7a2196fe`), 34 973 insertions and 16 705 deletions. The receipt-driven reviewer refuses the
whole candidate (`lens_context_budget_exceeded`): its reviewer evidence exceeds the native
context budget. Three independent START attempts on three different whole-branch targets hit
the same wall.

The mechanism that works is a **reduced `baseRef`**, not branch surgery. A review was started
for the last slice with an explicit base ref and the preflight passed, then the slice was
reviewed, approved and acknowledged:

```text
baseRef           e366678e7149ec87314068e1897367ec2f9f91bb
changed_files     2      original_changed_lines  212
risk_tier         medium selected_lenses         review-reliability
state             reviewing -> approved -> acknowledged (authority burned)
```

So each slice is reviewed with the tip of the previous slice as its base, with no branch
created and no history rewritten. The call shape is:

```text
gentle_review({"operation":"start","input":"{\"mode\":\"ordinary\",\"baseRef\":\"<previous tip>\",\"committedOnly\":true}"})
```

Slice 74 is the one already reviewed, approved and acknowledged.

## Slices

A slice never splits a commit. `lines` is insertions+deletions and is the size the reviewer
sees. `OVERSIZED` marks a single atomic commit that alone exceeds the budget and cannot be
reviewed without splitting that commit, which is a history rewrite.

| # | base-ref (previous tip) | tip | commits | lines | theme |
| --- | --- | --- | --- | --- | --- |
| 1 | `7a2196fe` | `04ecc78c` | 2 | 494 | chore(git): untrack local tooling state from .atl |
| 2 | `04ecc78c` | `59990b8e` | 1 | 8730 **OVERSIZED** | refactor(py): split the 4.3k-line PyO3 lib.rs into cohesiv |
| 3 | `59990b8e` | `aefafdc1` | 3 | 299 | docs(odd): record the 3D solid removal plan and its accept |
| 4 | `aefafdc1` | `1f843390` | 1 | 399 | test(ccx): make the CalculiX parity tests runnable and sto |
| 5 | `1f843390` | `473d6e89` | 2 | 416 | style(ccx): clear the lint violations in the CCX writer |
| 6 | `473d6e89` | `5baf8b81` | 1 | 2330 **OVERSIZED** | style: apply ruff format across the tree |
| 7 | `5baf8b81` | `9bb82441` | 1 | 35 | docs(odd): re-slice the solid removal along crate boundari |
| 8 | `9bb82441` | `7b295f3d` | 1 | 3457 **OVERSIZED** | refactor(rust): remove 3D solid element support |
| 9 | `7b295f3d` | `dd2e05aa` | 1 | 544 | refactor(py): remove the Python-side consumers of 3D solid |
| 10 | `dd2e05aa` | `9230ea2c` | 1 | 2292 **OVERSIZED** | refactor(py): remove the Python-side definitions of 3D sol |
| 11 | `9230ea2c` | `7e13d1e7` | 1 | 1829 **OVERSIZED** | test(ccx): delete the solid test paths and clean the docs |
| 12 | `7e13d1e7` | `5f751281` | 3 | 500 | fix(mesh): make the lazy gmsh binding static, and clear th |
| 13 | `5f751281` | `93f00119` | 4 | 279 | chore(py): clear the remaining src/ lint findings |
| 14 | `93f00119` | `13e9085a` | 2 | 489 | chore(tools): clear the last lint findings, leaving the no |
| 15 | `13e9085a` | `5b40ad59` | 1 | 204 | docs(odd): plan the test-suite physical correctness featur |
| 16 | `5b40ad59` | `4a46af29` | 2 | 487 | refactor(mitc4): delete the unwired drill-membrane operato |
| 17 | `4a46af29` | `e2526651` | 4 | 368 | test(parity): make the shell parity suite assert against t |
| 18 | `e2526651` | `ae343151` | 3 | 353 | test(mitc4): add the shell mass matrix invariants |
| 19 | `ae343151` | `5140d36b` | 3 | 457 | test(composites): replace the sign-blind B-coupling assert |
| 20 | `5140d36b` | `9cb49d14` | 2 | 254 | docs(odd): record the first verified code -> equation mapp |
| 21 | `9cb49d14` | `1b4c0684` | 4 | 454 | docs(odd): record the search for the most recent MITC3 wor |
| 22 | `1b4c0684` | `2fab0ba9` | 2 | 161 | docs(odd): plan the strain-smoothed MITC3+ implementation |
| 23 | `2fab0ba9` | `22b36517` | 1 | 571 | feat(mitc3): add the strain-smoothed membrane operator (S1 |
| 24 | `22b36517` | `80f50c32` | 1 | 645 | docs: add the shell element code-to-equation reference (U2 |
| 25 | `80f50c32` | `d003156c` | 5 | 460 | feat(mitc3): membrane-parameterised stiffness entry point |
| 26 | `d003156c` | `f8df4850` | 2 | 406 | style(mitc3): wrap two lines rustfmt wanted wrapped in the |
| 27 | `f8df4850` | `72e98f8b` | 3 | 478 | revert(mitc4): unwire the MITC4/D drill-membrane term, whi |
| 28 | `72e98f8b` | `2391beb3` | 2 | 246 | feat(mitc3): union-space force and tangent, and a global-D |
| 29 | `2391beb3` | `95cfd167` | 1 | 1040 | docs: add the materials and solver formulation references |
| 30 | `95cfd167` | `2eac7e66` | 5 | 315 | docs(rotor): drop the unverifiable ANSYS equation numbers |
| 31 | `2eac7e66` | `32df8fcf` | 1 | 630 | docs: add the validation matrix |
| 32 | `32df8fcf` | `688bb155` | 5 | 458 | docs: add Eqs. (10)-(16) to the MITC4+ extract |
| 33 | `688bb155` | `024baf71` | 3 | 470 | feat(mitc4): add the Winkler & Plakomytis ERC strain opera |
| 34 | `024baf71` | `45fc1314` | 1 | 479 | feat(mitc4): adopt the ERC as the element's drilling const |
| 35 | `45fc1314` | `fdfa531a` | 2 | 373 | test(ko2017): drop the xfail markers the ERC drilling fix |
| 36 | `fdfa531a` | `5ef5830d` | 1 | 192 | docs(mitc4plus): record what paper A actually defines, and |
| 37 | `5ef5830d` | `32299520` | 1 | 1906 **OVERSIZED** | docs(sdd): initialize the SDD change for the MITC4+/D elem |
| 38 | `32299520` | `6214bdb8` | 1 | 1709 **OVERSIZED** | docs(sdd): design, tasks and the 2025 formulation extract |
| 39 | `6214bdb8` | `e0107f26` | 1 | 245 | docs(mitc4plusd): transcribe paper B Eqs. (1)-(16) and set |
| 40 | `e0107f26` | `7f5370cd` | 1 | 658 | feat(mitc4_plusd): add the patch-test fixtures (SDD WU1) |
| 41 | `7f5370cd` | `daffa1d3` | 1 | 24 | docs(formulations): drop the overloaded "paper A/B" shorth |
| 42 | `daffa1d3` | `2f6cdf09` | 1 | 1158 | feat(mitc4_plusd): element core geometry, coefficients, di |
| 43 | `2f6cdf09` | `b50d1f10` | 1 | 1352 | feat(mitc4_plusd): the four B-operators, with the drill ve |
| 44 | `b50d1f10` | `ee132160` | 1 | 338 | feat(materials): additive uncorrected transverse-shear acc |
| 45 | `ee132160` | `d518cb02` | 1 | 1566 **OVERSIZED** | feat(mitc4_plusd): stiffness assembly, and fix a frame def |
| 46 | `d518cb02` | `8eec86c9` | 1 | 764 | docs(sdd): consolidate the apply-progress artifact |
| 47 | `8eec86c9` | `c5e800db` | 1 | 1091 | feat(mitc4_plusd): assembly-facing API, with the consisten |
| 48 | `c5e800db` | `8f525025` | 1 | 958 | test(mitc4_plusd): the 2017 paper's own basic tests, and t |
| 49 | `8f525025` | `bb0d671e` | 2 | 135 | docs(sdd): restate the zero-energy requirements to match p |
| 50 | `bb0d671e` | `82c7dc4f` | 1 | 508 | docs(sdd): second amendment to the zero-energy requirement |
| 51 | `82c7dc4f` | `acfca1de` | 4 | 336 | docs(sdd): Requirement 11's shearing scenario joins the in |
| 52 | `acfca1de` | `426249dc` | 2 | 420 | fix(mitc4): the production element's plane normal now uses |
| 53 | `426249dc` | `af8cf78d` | 1 | 1029 | test(mitc4_plusd): the 2025 paper's own basic tests, and T |
| 54 | `af8cf78d` | `3dfa7f6d` | 1 | 971 | test(mitc4_plusd): move the layout-bound Tier-2 tests onto |
| 55 | `3dfa7f6d` | `5e90b4bf` | 2 | 465 | fix(mitc4_plusd): the t-rule entered twice, making the E2 |
| 56 | `5e90b4bf` | `ea244c7c` | 2 | 387 | docs: handoff for the next session, and correct a wrong no |
| 57 | `ea244c7c` | `7ff506cd` | 1 | 7173 **OVERSIZED** | chore: checkpoint the uncommitted 2025-purity session (not |
| 58 | `7ff506cd` | `a0762753` | 1 | 356 | test(mitc4_plusd): reference-free diagnostics that localis |
| 59 | `a0762753` | `1dad1ebc` | 3 | 361 | docs(mitc4_plusd): the round-off bound is 1.1e-5, not 5e-6 |
| 60 | `1dad1ebc` | `923260e1` | 11 | 489 | test(modal): a non-vacuity guard for the Rust-vs-Python mo |
| 61 | `923260e1` | `70f1e485` | 4 | 249 | docs(odd): save the modal participation-factor verificatio |
| 62 | `70f1e485` | `f555d13a` | 2 | 459 | docs(odd): iteration 24 -- the analytic B/N route specifie |
| 63 | `f555d13a` | `71838ce9` | 5 | 404 | docs(odd): wiring attempt #1 was structurally wrong (116%) |
| 64 | `71838ce9` | `bbc39888` | 5 | 339 | test(mitc4_plusd): the square fixture refutes the metric r |
| 65 | `bbc39888` | `56bfd95b` | 1 | 1204 | test(suite): zero ignored and zero failures, with the limi |
| 66 | `56bfd95b` | `0c0d91df` | 1 | 26 | docs(openspec): close the bookkeeping -- 10.6 done, 13.1 b |
| 67 | `0c0d91df` | `b479c38f` | 1 | 23084 **OVERSIZED** | refactor(mitc4): retire the hybrid and take the name mitc4 |
| 68 | `b479c38f` | `b097549c` | 1 | 369 | docs(formulations): audit the element documentation agains |
| 69 | `b097549c` | `133f40cb` | 4 | 492 | docs(openspec): close 12.3, 12.4 and record the 13.1 close |
| 70 | `133f40cb` | `2fd847dc` | 1 | 569 | chore(openspec): archive mitc4plusd-faithful and apply its |
| 71 | `2fd847dc` | `15215c03` | 3 | 366 | docs(odd): reconcile the formulation-documentation task wi |
| 72 | `15215c03` | `a13c865a` | 1 | 373 | test(composite): parity of five layups against CalculiX S8 |
| 73 | `a13c865a` | `98ba3094` | 3 | 480 | docs(validation-matrix): add the composite layup and IEA b |
| 74 | `98ba3094` | `3fefba68` | 6 | 497 | docs(validation-matrix): add the blade static rows and the |

## What still blocks

10 of the 74 slices are oversized single commits:

- slice 2, tip `59990b8e`, 8730 lines — refactor(py): split the 4.3k-line PyO3 lib.rs into cohesive modules
- slice 6, tip `5baf8b81`, 2330 lines — style: apply ruff format across the tree
- slice 8, tip `7b295f3d`, 3457 lines — refactor(rust): remove 3D solid element support
- slice 10, tip `9230ea2c`, 2292 lines — refactor(py): remove the Python-side definitions of 3D solid support
- slice 11, tip `7e13d1e7`, 1829 lines — test(ccx): delete the solid test paths and clean the docs
- slice 37, tip `32299520`, 1906 lines — docs(sdd): initialize the SDD change for the MITC4+/D element
- slice 38, tip `6214bdb8`, 1709 lines — docs(sdd): design, tasks and the 2025 formulation extract
- slice 45, tip `d518cb02`, 1566 lines — feat(mitc4_plusd): stiffness assembly, and fix a frame defect it exposed (SDD WU4)
- slice 57, tip `7ff506cd`, 7173 lines — chore: checkpoint the uncommitted 2025-purity session (not delivery-ready)
- slice 67, tip `b479c38f`, 23084 lines — refactor(mitc4): retire the hybrid and take the name mitc4 (SDD 11.1-11.5)

The largest, `b479c38` (slice 67, 23 084 lines), is the retirement of the hybrid MITC4+ and
the move of the faithful element into `mitc4.rs`: a file move plus a rewrite, so its diff is
the whole element. Reviewing it below the budget means splitting that commit itself, which is
a history rewrite of the branch, not a slice of it.

## Recommendation

1. Review the slices under the budget in order, each with `baseRef` = the previous tip. No
   branch is created and the working tree is untouched; the current branch stays the source
   of truth.
2. For the oversized commits, decide per commit: split it into smaller commits on a review
   branch (keeps the current branch intact, but the split is real work), or record it as a
   grandfathered baseline by explicit maintainer exception.
3. Do not re-run START on the whole-branch candidate: it is documented to fail, three times
   over, and the provider marks it `not_replayable` / `stop`.
