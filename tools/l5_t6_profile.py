"""T15 helper: time the Python pieces of the rotor step callback (no production change).

Usage: scripts/aeroenv.sh python tools/l5_t6_profile.py <solid.yaml> --workdir <dir>
Monkeypatches timers around the methods `_step_cb` calls and prints a summary at exit.
"""
import atexit
import sys
import time
from collections import defaultdict

from aeroelast.solvers.fsi import linear_dynamic, rotor

T = defaultdict(lambda: [0.0, 0])


def wrap(cls, name):
    f = getattr(cls, name, None)
    if f is None:
        return

    def g(*a, **k):
        t0 = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            d = T[f"{cls.__name__}.{name}"]
            d[0] += time.perf_counter() - t0
            d[1] += 1

    setattr(cls, name, g)


for cls in (rotor.LinearDynamicFSIRotorCorotationalSolver, linear_dynamic.LinearDynamicFSISolver):
    for n in ("_compute_stress_fields", "_log_rotor_performance", "_write_restart_state",
              "_log_structural_report", "_log_probe_data", "_handle_checkpoint"):
        wrap(cls, n)


@atexit.register
def _report():
    with open("t6_profile_summary.txt", "w") as fh:
        for k, (s, n) in sorted(T.items(), key=lambda kv: -kv[1][0]):
            fh.write(f"{k}: total {s:.2f} s over {n} calls = {s / max(n, 1) * 1e3:.0f} ms/call\n")


from aeroelast.cli.aeroelast import main  # noqa: E402

sys.exit(main())
