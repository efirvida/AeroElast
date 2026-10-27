# T21 finding: CCX `tip_axial_tension` work metric does NOT meet 5 %

Revision: worktree input-fix-CCXW, branch fix/ccx-parity-work, base ad5f3ac + uncommitted diff to tests/test_blade_ccx_parity.py.
Metric: W = 1/2 sum f.u (tip_axial_tension, 1e6 N axial), AE vs CalculiX S8R reference. Bound: 5 % (not widened).

| tree | W_AE | W_CCX | work gap | old node-mean gap | verdict |
|---|---|---|---|---|---|
| T7e  | 3.256675e+05 | 3.600873e+05 |  9.56 % | 4.64 % | FAIL |
| T8w  | 2.791395e+05 | 3.123260e+05 | 10.63 % | 4.63 % | FAIL |
| T14m (final-tree candidate, "combination") | 4.045477e+05 | 5.090555e+05 | 20.53 % | 7.07 % | FAIL |

Jobs: 11623281 (T7e), 11623282 (T8w), 11623283 (T14m); logs work-w-<tree>-tip_axial_tension.log.
11623284 (full CCX parity) hit TIMEOUT at 1:30; no per-case results. Not resubmitted: the gating case already fails, other load cases cannot change the verdict.
Unit tests (test_parity_metric_ignores_node_count, test_load_work_value_and_missing_node): 2 passed; ruff clean.

Observation: the work gap is larger than the node-mean gap on every tree (MITC4 is softer in W by 10-20 %), so the work metric does not rescue the 5 % band; it exposes it. W is not independent of tip compliance, and T14m differs from CCX the most.

Verdict: FAIL on all three trees and on the combination.
Recommendation: record as finding; do not commit/merge fix/ccx-parity-work, change nothing else. Next step is to understand the 20 % W gap on T14m (T14m mesh differs from T7e/T8w), not to alter the tolerance.
