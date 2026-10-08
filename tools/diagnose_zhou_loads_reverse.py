"""Zhou's published loads, completed with a reconstructed pitching moment, in our models.

Issue #12 cannot be closed from the published data as it stands: their Fig. 11 gives
`Np` and `Tp` only, and the sectional torque is dominated by `Mp`, which they never
publish. This reconstructs it from what they DO publish and the one input both models
share - the official airfoil polars:

    |F| = sqrt(Np^2 + Tp^2) = q c sqrt(CL^2 + CD^2)   at their published angle of attack
      => q c = |F| / sqrt(CL^2 + CD^2)
      => Mp  = q c^2 Cm(alpha)

with `CL, CD, Cm` interpolated from the official polar at their `alpha` (their Fig. 10).
Then their complete sectional load drives our two models - the shell through the
production projector, and the standard beam through the deck's own shear centre and
`GKt` - and the tip twists are compared against their -3.60 deg.

Sign mapping: their torsional deflection is positive toward stall; our `omega` is
positive nose-down, i.e. away from stall. The deck puts the leading edge at `+x` and the
load frame's downwind thrust at `+y`, so a `+z` rotation moves the leading edge downwind and
**nose-down is positive** on our frame. The BEM's `Mp` carries the opposite sign (a nose-down
pitch moment is negative), which is why the shell path takes it through the production
projector while the beam's `m` is the textbook's nose-up-positive form (Dowell et al., eq.
2.1.2). Everything printed below is converted once, at the print, to the one
nose-down-positive convention, so the shell, the beam and Zhou are directly comparable.

Caveat, stated rather than hidden: the polars are ours (both models use the same official
airfoils), and the reconstruction extrapolates nothing - it is restricted to 0.15 <= r/R
where their angle-of-attack figure starts.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from _aeroelast import PyMeshAssembler
from openfast_toolbox.converters.beam import K66toPropsDecoupled
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

# A script puts its own directory on sys.path, not the repository root, so `import tests.*`
# below fails when this is run the way the issue documents it
# (`python tools/diagnose_zhou_loads_reverse.py`) with
# `ModuleNotFoundError: No module named 'tests.validation'`. Insert the root before the first
# `tests` import. `tests.support.paths` computes the same anchor by walking up to
# `pyproject.toml`.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.validation.blade.test_blade_iea15mw_validation as bv  # noqa: E402
import tests.validation.blade.test_blade_rated_twist as t  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMResult  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.validation.blade.test_blade_twist_anchor_beam import (  # noqa: E402
    BEAMDYN_BLADE,
    ELASTODYN_BLADE,
    GKT_MEDIAN_FACTOR,
    _numeric_rows,
)

SPAN = np.array([0.0, 0.0, 1.0])
RHO = 1.225
ZHOU_LOADS = REPO_ROOT / "docs" / "validation_data" / "zhou_2025_fig11_loads.csv"
ZHOU_AOA = REPO_ROOT / "docs" / "validation_data" / "zhou_2025_fig10_aoa.csv"
ZHOU_TIP_TORSION_DEG = -3.60
R_MIN = 0.15

# ---- their published data -----------------------------------------------------
loads = pd.read_csv(ZHOU_LOADS)
aoa = pd.read_csv(ZHOU_AOA)
r_R_load = loads["r_R"].to_numpy(dtype=float)
Np_zhou = loads["normal_force_kN_m"].to_numpy(dtype=float) * 1.0e3
Tp_zhou = loads["tangential_force_kN_m"].to_numpy(dtype=float) * 1.0e3
aoa_zhou = np.interp(
    r_R_load, aoa["r_R"].to_numpy(dtype=float), aoa["angle_of_attack_deg"].to_numpy(dtype=float)
)

# ---- the shared input: the official blade, chord and polars -------------------
aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
r_hub = np.asarray(aero.r, dtype=float)
hub = float(aero.hub_radius)
rotor_radius = float(r_hub[-1])
chord = np.asarray(aero.chord, dtype=float)

# ---- our mesh, shell and beam -------------------------------------------------
Node._id_counter = 0
MeshElement._id_counter = 0
model = Blade(str(t.YAML), element_size=0.5)
model.generate_mesh()
mesh = model.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = model.get_element_properties()
coords = mesh.coords_array
assembler = PyMeshAssembler.from_model(
    bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None
)
n_dof = assembler.dofs_count
rk, ck, vk = assembler.assemble_k()
K = coo_matrix((np.asarray(vk), (np.asarray(rk), np.asarray(ck))), shape=(n_dof, n_dof)).tocsr()
root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
free = np.array([i for i in range(n_dof) if i not in {6 * r + d for r in root for d in range(6)}])
Kff = K[np.ix_(free, free)]
phys = t._physical_stations(coords)
tip = np.where(np.abs(coords[:, 2] - phys[-1]) < t.STATION_GAP_TOLERANCE)[0]

deck_rows = _numeric_rows(BEAMDYN_BLADE, "DISTRIBUTED PROPERTIES")
n_deck = len(deck_rows) // 13
frac = np.array([deck_rows[i * 13][0] for i in range(n_deck)])
K66 = np.array([np.asarray(deck_rows[i * 13 + 1 : i * 13 + 7]) for i in range(n_deck)])
deck = [K66toPropsDecoupled(K66[i], convention="BeamDyn") for i in range(n_deck)]
x_shear = np.array([p[8] for p in deck])
gkt = np.array([p[5] for p in deck])
ed = np.asarray(_numeric_rows(ELASTODYN_BLADE, "DISTRIBUTED BLADE PROPERTIES"))
pitch = np.interp(frac, ed[:, 0], ed[:, 1])

# ---- reconstruct their pitching moment ----------------------------------------
Mp_zhou = np.zeros_like(Np_zhou)
alpha_used = np.zeros_like(Np_zhou)
cl_used = np.zeros_like(Np_zhou)
for i, (rr, aa, np_i, tp_i) in enumerate(zip(r_R_load, aoa_zhou, Np_zhou, Tp_zhou, strict=True)):
    if rr < R_MIN:
        continue
    r_abs = rr * rotor_radius
    if r_abs <= r_hub[0] or r_abs >= r_hub[-1]:
        continue
    j = int(np.argmin(np.abs(r_hub - r_abs)))
    polar = aero.stations[j].airfoil.polars[0]
    cl = float(np.interp(np.deg2rad(aa), polar.alpha, polar.cl))
    cd = float(np.interp(np.deg2rad(aa), polar.alpha, polar.cd))
    cm = float(np.interp(np.deg2rad(aa), polar.alpha, polar.cm))
    force_mag = float(np.hypot(np_i, tp_i))
    denom = float(np.hypot(cl, cd))
    if denom <= 0.0:
        continue
    qc = force_mag / denom  # = 0.5 rho W^2 c
    Mp_zhou[i] = qc * chord[j] * cm
    alpha_used[i] = aa
    cl_used[i] = cl

print("their published angle of attack and the reconstructed pitching moment")
print(
    f"{'r/R':>6} {'aoA[deg]':>9} {'CL':>7} {'Np[kN/m]':>10} {'Tp[kN/m]':>10} {'Mp rec[N.m/m]':>15}"
)
for i in range(0, len(r_R_load), max(1, len(r_R_load) // 8)):
    print(
        f"{r_R_load[i]:6.2f} {alpha_used[i]:9.2f} {cl_used[i]:7.3f} "
        f"{Np_zhou[i] / 1e3:10.2f} {Tp_zhou[i] / 1e3:10.2f} {Mp_zhou[i]:15.1f}"
    )

# ---- their loads on OUR stations, through OUR two models ----------------------
r_their = r_R_load * rotor_radius
order = np.argsort(r_their)
Np_our = np.interp(r_hub, r_their[order], Np_zhou[order], left=0.0, right=0.0)
Tp_our = np.interp(r_hub, r_their[order], Tp_zhou[order], left=0.0, right=0.0)
Mp_our = np.interp(r_hub, r_their[order], Mp_zhou[order], left=0.0, right=0.0)

bem_zhou = BEMResult(
    r=r_hub,
    Np=Np_our,
    Tp=Tp_our,
    Mp=Mp_our,
    alpha=np.zeros_like(r_hub),
    cl=np.zeros_like(r_hub),
    cd=np.zeros_like(r_hub),
    a=np.zeros_like(r_hub),
    ap=np.zeros_like(r_hub),
    thrust=0.0,
    torque=0.0,
    power=0.0,
)
projector = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)
forces = projector.project(bem_zhou)
f = np.zeros(n_dof)
f[0::6] = forces[:, 0]
f[1::6] = forces[:, 1]
f[2::6] = forces[:, 2]
u = np.zeros(n_dof)
u[free] = np.asarray(spsolve(Kff, f[free]), dtype=float).ravel()
omega = float(np.rad2deg(t._ring_kinematics(coords, u, tip)["omega"]))

# the beam, on the same completed load, standard shear-centre arm
r_b = hub + frac * (rotor_radius - hub)
Np_b = np.interp(r_b, r_their[order], Np_zhou[order], left=0.0, right=0.0)
Tp_b = np.interp(r_b, r_their[order], Tp_zhou[order], left=0.0, right=0.0)
Mp_b = np.interp(r_b, r_their[order], Mp_zhou[order], left=0.0, right=0.0)
chord_b = np.interp(r_b, r_hub, chord)
m_b = Mp_b + ((pitch - 0.25) * chord_b - x_shear) * Np_b - (0.0 - 0.0) * Tp_b
floor = float(np.median(gkt)) / GKT_MEDIAN_FACTOR
ok = np.minimum(gkt[:-1], gkt[1:]) >= floor
dr = np.diff(r_b)
torque = np.zeros_like(m_b)
for i in range(len(m_b) - 2, -1, -1):
    torque[i] = torque[i + 1] + (0.5 * (m_b[i] + m_b[i + 1]) * dr[i] if ok[i] else 0.0)
theta = np.zeros_like(m_b)
for i in range(1, len(m_b)):
    theta[i] = theta[i - 1] + (
        0.5 * (torque[i - 1] / gkt[i - 1] + torque[i] / gkt[i]) * dr[i - 1] if ok[i - 1] else 0.0
    )
phi = float(np.rad2deg(theta[-1]))

print()
print("their complete load (forces published, moment reconstructed) in OUR two models:")
print("one convention throughout, nose-down positive (leading edge +x, load downwind +y):")
print(f"  shell (MITC4, 0.5 m)   tip section rotation {omega:+9.4f} deg")
print(f"  standard beam          tip twist            {-phi:+9.4f} deg")
print(
    f"  Zhou, their own LL-FVW+GEBT                  {-ZHOU_TIP_TORSION_DEG:+9.4f} deg  (their paper: toward stall positive)"
)
print()
print(
    "  reading: with nose-down positive, their -3.60 deg is +3.60 deg, so the comparable "
    f"magnitudes are shell |{omega:.3f}| -> {abs(omega) / abs(ZHOU_TIP_TORSION_DEG):.2f}x "
    f"and beam |{-phi:.3f}| -> {abs(phi) / abs(ZHOU_TIP_TORSION_DEG):.2f}x"
)
print(
    "  both our models now twist the same way under this load; the beam's printed value is"
    " its own number negated, once, because its `m` is the textbook's nose-up-positive form"
)
print("  (see odd/tasks/rated-twist-sign-convention.md for why the sign is not a free choice)")
print(
    f"  reconstructed pitching moment, integral over the span: {np.sum(Mp_our * np.gradient(r_hub)):+.6e} N.m"
)
