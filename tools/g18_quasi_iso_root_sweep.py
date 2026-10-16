"""G18 quasi_iso transverse bending: is the ~1.0-1.2 % offset the clamped-root restraint?

Predeclared criterion (written before the first run). The hand-CLT strip closed form assumes a
free-width strip; a clamped root restrains the width-wise curvature and the twist-bending coupling
(D16, D26, non-zero for the quasi-isotropic layup) over a length that scales with the width, so
the relative effect scales with B/L. AeroElast alone (the CalculiX comparison is non-gating) is run
at B = 0.1, 0.05, 0.025, 0.0125 m with L = 1 m, ny = 80 and nx in {4, 8} (the second mesh checks
that the gap is not a discretization effect). The offset is ATTRIBUTED to root restraint when
(a) the gap at B = 0.1 reproduces the stored -1 to -1.2 % band, (b) |gap| decreases monotonically
as B/L decreases, (c) |gap| at B = 0.0125 is <= 0.3 % (the level of the other layups), and (d) the
gap changes < 0.3 percentage points between nx = 4 and 8 at every B, and (e) sym_0_90s closes
at <= 0.3 % on both meshes at every width. This control gates attribution. Run through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import sys

import numpy as np

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_layup_parity as m  # noqa: E402
from test_composite_beam_parity import _hand_clt_abd  # noqa: E402


def analytic_bending(angles):
    A, Bm, D = _hand_clt_abd(m.E1, m.E2, m.G12, m.NU12, list(angles), m.THICKNESS / len(angles))
    inv = np.linalg.inv(np.block([[A, Bm], [Bm, D]]))
    return inv[4, 4] * 100.0 * m.L**3 / (3.0 * m.B)


dof, load, mdof = m.CASES["bending"]
gaps = {}
for name in ("quasi_iso", "sym_0_90s"):
    angles = m.LAYUPS[name]
    lam = m._laminate(angles)
    for width in (0.1, 0.05, 0.025, 0.0125):
        m.B = width
        ref = analytic_bending(angles)
        row = []
        for nx in (4, 8):
            mesh = m._build_mesh(nx, 80)
            face = [mesh.node_id_to_index[n.id] for n in mesh.get_node_set("free_face").nodes.values()]
            u = m._aero_solution(mesh, lam, dof, load)
            w = float(np.mean([u[6 * k + mdof] for k in face]))
            row.append(100 * (w / ref - 1))
        gaps[(name, width)] = row
        print(f"{name:10s} B={width:<7} B/L={width/m.L:<7} gap nx=4 {row[0]:+7.3f} %  nx=8 {row[1]:+7.3f} %", flush=True)

q = [gaps[("quasi_iso", w)] for w in (0.1, 0.05, 0.025, 0.0125)]
a = -1.2 <= q[0][1] <= -1.0
b = all(abs(q[i + 1][1]) < abs(q[i][1]) for i in range(3))
c = abs(q[3][1]) <= 0.3
d = all(abs(r[0] - r[1]) < 0.3 for r in q)
e = all(abs(gap) <= 0.3 for (name, _), row in gaps.items() if name == "sym_0_90s" for gap in row)
print(f"criteria: (a) B=0.1 in band {a}  (b) monotone {b}  (c) B=0.0125 <= 0.3 % {c}  (d) mesh {d}  (e) control <= 0.3 % {e}")
print("  => " + ("ATTRIBUTED" if all((a, b, c, d, e)) else "UNATTRIBUTED") + " under the predeclared criterion")

# The test's own mesh is 4x10: show how much of its stored gap is length-wise discretization.
m.B = 0.1
lam, ref = m._laminate(m.LAYUPS["quasi_iso"]), analytic_bending(m.LAYUPS["quasi_iso"])
for ny in (10, 20, 40, 80):
    mesh = m._build_mesh(4, ny)
    face = [mesh.node_id_to_index[n.id] for n in mesh.get_node_set("free_face").nodes.values()]
    u = m._aero_solution(mesh, lam, dof, load)
    print(f"quasi_iso B=0.1 nx=4 ny={ny:<3} gap {100 * (np.mean([u[6 * k + mdof] for k in face]) / ref - 1):+.3f} %")
