"""Invert the twist: what torsional moment does the reference need, and do we have it?

Issue #12's first work unit, run with the corrected sign chain: instead of asking
"is our twist short", it asks what the reference implies about the load.

The spanwise sectional torsional moment is taken about the section's **shear
centre**, in the deck's own frame:

    m(z) = Mp + (x_ac - xS) Np - (y_ac - yS) Tp,   x_ac = (pitch - 0.25) chord

with `pitch` the ElastoDyn pitch axis, `xS, yS, GKt` from the BeamDyn blade deck
through the same `K66toPropsDecoupled` route and the same thin-tip guard
`test_blade_twist_anchor_beam.py` uses, so this tool and that test cannot disagree
about the deck or the integration. The tip twist it returns is the cross-check
against that test's own number.

Measured: the physical arm (aerodynamic centre relative to the shear centre) is
about 0.06 c, where the projector transfers the moment from a reference 0.22 c
behind the aerodynamic centre - the arithmetic-mean node centroid of the mesh.
"""

import numpy as np

import tests.validation.blade.test_blade_rated_twist as t
from openfast_toolbox.converters.beam import K66toPropsDecoupled
from tests.support.openfast_bem import build_blade_aero_from_aerodyn
from tests.validation.blade.test_blade_twist_anchor_beam import (
    BEAMDYN_BLADE,
    ELASTODYN_BLADE,
    GKT_MEDIAN_FACTOR,
    _numeric_rows,
)
from aeroelast.solvers.bem.engine import BEMSolver

ZHOU_TIP_TORSION_DEG = -3.60
POLLARS_AC = 0.25

rows = _numeric_rows(BEAMDYN_BLADE, "DISTRIBUTED PROPERTIES")
n_st = len(rows) // 13
frac = np.array([rows[i * 13][0] for i in range(n_st)])
K = np.array([np.asarray(rows[i * 13 + 1 : i * 13 + 7]) for i in range(n_st)])
deck = [K66toPropsDecoupled(K[i], convention="BeamDyn") for i in range(n_st)]
x_shear = np.array([p[8] for p in deck])
y_shear = np.array([p[9] for p in deck])
gkt = np.array([p[5] for p in deck])

elastodyn = np.asarray(_numeric_rows(ELASTODYN_BLADE, "DISTRIBUTED BLADE PROPERTIES"))
pitch_frac, pitch_axis = elastodyn[:, 0], elastodyn[:, 1]

aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
bem = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
    t.V_RATED, t.RPM_RATED, t.PITCH_RATED
)
if bem.Mp is None:
    raise SystemExit("no pitching moment")

r_hub = np.asarray(aero.r, dtype=float)
r = float(aero.hub_radius) + frac * float(aero.blade_length)
Np = np.interp(r, r_hub, bem.Np)
Tp = np.interp(r, r_hub, bem.Tp)
Mp = np.interp(r, r_hub, bem.Mp)
chord = np.interp(r, r_hub, aero.chord)
pitch = np.interp(frac, pitch_frac, pitch_axis)

x_ac = (pitch - POLLARS_AC) * chord  # AC offset from the pitch axis, y_ac = 0
arm = (x_ac - x_shear) * Np - (0.0 - y_shear) * Tp
m = Mp + arm

# The anchor test's thin-tip guard: an interval whose GKt collapses is not integrated.
floor = float(np.median(gkt)) / GKT_MEDIAN_FACTOR
interval_ok = np.minimum(gkt[:-1], gkt[1:]) >= floor
dr = np.diff(r)

torque = np.zeros_like(m)  # T(z) = int_z^R m ds, root-fixed
for i in range(len(m) - 2, -1, -1):
    torque[i] = torque[i + 1] + (0.5 * (m[i] + m[i + 1]) * dr[i] if interval_ok[i] else 0.0)
twist = np.zeros_like(m)  # theta(z) = int_0^z T/GJ ds
for i in range(1, len(m)):
    step = (
        0.5 * (torque[i - 1] / gkt[i - 1] + torque[i] / gkt[i]) * dr[i - 1]
        if interval_ok[i - 1]
        else 0.0
    )
    twist[i] = twist[i - 1] + step

print(f"anchor stations: {n_st}   blade length {r[-1] - r[0]:.3f} m")
print(f"physical arm x_ac - xS range          = {(x_ac - x_shear).min():+.4f} .. {(x_ac - x_shear).max():+.4f} m")
print(f"  as a chord fraction                 = {((x_ac - x_shear) / chord).min():+.4f} .. {((x_ac - x_shear) / chord).max():+.4f} c")
print()
print(f"integral of m(z) over the span        = {np.sum(0.5 * (m[:-1] + m[1:]) * dr):+.6e} N.m")
print(f"  of which the polars' Mp             = {np.sum(0.5 * (Mp[:-1] + Mp[1:]) * dr):+.6e} N.m")
print(f"  of which the shear-centre arm       = {np.sum(0.5 * (arm[:-1] + arm[1:]) * dr):+.6e} N.m")
print()
print(f"beam tip twist with this m(z)         = {np.rad2deg(twist[-1]):+.4f} deg   (anchor test: -1.7765)")
print(f"Zhou's tip torsion                    = {ZHOU_TIP_TORSION_DEG:+.2f} deg")
print(f"scale our m(z) would need             = {ZHOU_TIP_TORSION_DEG / np.rad2deg(twist[-1]):+.3f}x")
print()
print(f"{'r[m]':>8} {'m[N.m/m]':>12} {'Mp':>12} {'arm':>12} {'x_ac-xS[m]':>11} {'GKt[N.m2]':>12}")
for i in range(0, n_st, max(1, n_st // 10)):
    print(
        f"{r[i]:8.2f} {m[i]:12.4e} {Mp[i]:12.4e} {arm[i]:12.4e} "
        f"{x_ac[i] - x_shear[i]:11.4f} {gkt[i]:12.4e}"
    )
