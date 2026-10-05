"""What section torque are we actually applying at rated, and which way does it point?

The multi-cell realisation is verified correct (skin walls carry q_i, shared webs
q_i - q_j, total moment exact, net force zero) yet it gives a nose-up-ish tip
rotation (+0.0773 deg) where the minimum-norm field gives nose-down (-0.8696 deg)
for the *same* requested torque. Two fields with the same resultant moment cannot
differ in the sign of a rigid rotation, so the number to check is the requested
torque itself:

    torsion_k = M_strip_k . span_hat
              = Mp_k * dr_k  -  cross(ac_offset_k, F_strip_k) . span_hat

i.e. the aerodynamic pitching moment versus the lever-arm transfer, and their
signs. AeroDyn's rated pitching moment is nose-down (negative), so a positive
total would mean the transfer term is overwhelming it with the wrong sign.
"""

import numpy as np

import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.engine import BEMSolver
from aeroelast.solvers.bem.force_projection import ForceProjector
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])

blade = Blade(str(t.YAML), element_size=1.0)
blade.generate_mesh()
mesh = blade.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = blade.get_element_properties()
aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
proj = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)

solver = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0)
bem = solver.compute(t.V_RATED, t.RPM_RATED, t.PITCH_RATED)
Mp = bem.Mp
if Mp is None:
    raise SystemExit("the rated BEM result carries no pitching moment")

m_mp = []
m_transfer = []
for k, strip in enumerate(proj._strips):
    F_strip = (
        float(bem.Np[k]) * strip.dr * proj._strip_normal_dirs[k]
        + float(bem.Tp[k]) * strip.dr * proj._strip_chord_dirs[k]
    )
    m_mp.append(float(Mp[k]) * strip.dr)
    m_transfer.append(-float(np.cross(proj._strip_ac_offsets[k], F_strip) @ SPAN))

m_mp = np.asarray(m_mp)
m_transfer = np.asarray(m_transfer)
total = m_mp + m_transfer
r = np.asarray(aero.r)
i_mp = float(np.sum(0.5 * (Mp[:-1] + Mp[1:]) * np.diff(r)))

print(f"rated: V={t.V_RATED} m/s, {t.RPM_RATED} rpm")
print(
    f"BEM pitching moment  int(Mp dr)      = {i_mp:+.6e} N.m   (per blade; nose-down if negative)"
)
print(f"Mp  * dr, summed over strips         = {m_mp.sum():+.6e} N.m")
print(f"transfer -cross(ac_offset,F).span    = {m_transfer.sum():+.6e} N.m")
print(f"TOTAL applied section torque         = {total.sum():+.6e} N.m")
print(f"  -> sign: {'NOSE-UP (positive)' if total.sum() > 0 else 'NOSE-DOWN (negative)'}")
print()
print("where is the node-mean (the moment reference) along the chord?")
print(f"  {'k':>3} {'c[m]':>7} {'x_mean/c':>9} {'x_AC/c':>7} {'transfer sign':>14}")
for k in range(0, len(proj._strips), max(1, len(proj._strips) // 8)):
    strip = proj._strips[k]
    if len(strip.node_indices) < 3:
        continue
    ring = next(
        (g for g in proj._strip_ring_groups[k] if len(g) >= 3),
        None,
    )
    if ring is None:
        continue
    pts = strip.centroid + strip.offsets[ring]
    ch = proj._strip_chord_dirs[k]
    le_i, te_i = ForceProjector._section_ends(pts, ch, SPAN)
    chord = float(np.linalg.norm(pts[le_i] - pts[te_i]))
    le = pts[le_i]
    x_mean = float((strip.centroid - le) @ ch) / chord
    ac_frac = float(aero.stations[k].airfoil.aerodynamic_center)
    t_sign = "NOSE-UP" if m_transfer[k] > 0 else "nose-down"
    print(f"  {k:3d} {chord:7.3f} {x_mean:9.3f} {ac_frac:7.3f} {t_sign:>14}")
print()
print("per-strip breakdown (r, Mp*dr, transfer, total):")
for k in range(0, len(total), max(1, len(total) // 8)):
    print(
        f"  strip {k:2d}  r={r[k]:7.2f}  Mp*dr={m_mp[k]:+12.4e}  "
        f"transfer={m_transfer[k]:+12.4e}  total={total[k]:+12.4e}"
    )
