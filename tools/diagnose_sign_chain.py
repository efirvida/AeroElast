"""Arbitrate the sign chain against the deck's own aerofoil geometry.

Three conventions meet in the applied torsional load and none was pinned to an
external reference:

1. which end of a mesh ring is the leading edge (`_section_ends` guesses from the
   blunt end);
2. the sense of the chord axis `_strip_chord_dirs`, resolved from the configured
   tangential direction, which then decides the axis `Mp` is applied on as
   `Mp * dr * span_hat`;
3. the sign convention of the section-rotation estimator that reports "nose-down".

The deck settles all three. The canonical WindIO blade (`tests/IEA-15-240-RWT.yaml`)
carries every aerofoil's coordinates with x = 0 at the leading edge, and the mesh was
built from exactly those tables, so matching a ring to its aerofoil identifies the
leading edge without the projector's opinion. This prints, per station:

  * the aerofoil name the station maps to and the match residual for each of the
    four possible alignments, so the winner is not a guess;
  * whether the deck's leading-to-trailing-edge direction runs along
    `_strip_chord_dirs` (+1) or against it (-1) - the factor `Mp` needs;
  * which x side the deck's leading edge sits on, and the estimator's sign for a
    rigid rotation about +z, hence whether the module's "omega < 0 is nose-down"
    holds on this frame.
"""

import numpy as np
import yaml

import tests.validation.blade.test_blade_rated_twist as t
from aeroelast.models.blade.model import Blade
from aeroelast.solvers.bem.force_projection import ForceProjector
from tests.support.openfast_bem import build_blade_aero_from_aerodyn

SPAN = np.array([0.0, 0.0, 1.0])
YAML = t.YAML

deck = yaml.safe_load(YAML.read_text())
airfoils = {
    a["name"]: (np.asarray(a["coordinates"]["x"], float), np.asarray(a["coordinates"]["y"], float))
    for a in deck["airfoils"]
}
position = deck["components"]["blade"]["outer_shape_bem"]["airfoil_position"]
labels = position["labels"]
grid = np.asarray(position["grid"], float)

blade = Blade(str(YAML), element_size=1.0)
blade.generate_mesh()
mesh = blade.mesh
if mesh is None:
    raise SystemExit("no mesh")
props = blade.get_element_properties()
coords = mesh.coords_array
aero = build_blade_aero_from_aerodyn(t.AD_PRIMARY)
proj = ForceProjector(mesh, aero, span_direction=SPAN, element_properties=props)

hub = float(aero.hub_radius)
length = float(aero.r[-1] - aero.r[0])


def deck_airfoil(fraction: float) -> str:
    i = int(np.argmin(np.abs(grid - fraction)))
    return labels[i]


print(
    f"{'k':>3} {'af':>16} {'resid/c':>8} {'LE . c^':>9} {'Mp factor':>10} {'LE at x':>8} {'_section_ends':>14} {'agree':>6}"
)
for k in range(2, len(proj._strips), max(1, len(proj._strips) // 8)):
    strip = proj._strips[k]
    ring = next((g for g in proj._strip_ring_groups[k] if len(g) >= 3), None)
    if ring is None:
        continue
    pts = strip.centroid + strip.offsets[ring]
    c_hat = proj._strip_chord_dirs[k]
    n_hat = proj._strip_normal_dirs[k]
    in_plane = np.cross(SPAN, c_hat)
    in_plane /= np.linalg.norm(in_plane)

    fraction = float((aero.r[k] - hub) / length)
    name = deck_airfoil(fraction)
    ax, ay = airfoils[name]
    x_proj = pts @ c_hat
    y_proj = pts @ in_plane
    chord = float(x_proj.max() - x_proj.min())

    # align the normalised aerofoil to the ring: 4 flips, start at the ring's low end
    best = None
    for flip_x in (1.0, -1.0):
        for flip_y in (1.0, -1.0):
            qx = x_proj.min() + chord * (1.0 - ax if flip_x < 0 else ax)
            qy = y_proj.mean() + chord * flip_y * ay
            resid = []
            for xi, yi in zip(x_proj, y_proj, strict=True):
                d = np.hypot(qx - xi, qy - yi)
                resid.append(float(d.min()))
            score = float(np.mean(resid)) / chord
            cand = (score, flip_x, flip_y)
            if best is None or score < best[0]:
                best = cand
    score, flip_x, flip_y = best

    # deck leading edge direction in the ring's own chord frame
    le_at_low = flip_x > 0  # x=0 (the leading edge) sits at the low projected end
    # direction from leading to trailing edge along c_hat
    le_to_te = +1.0 if le_at_low else -1.0
    mp_factor = +1.0 if le_to_te > 0 else -1.0

    # which x side is the deck leading edge on
    le_index = int(np.argmin(np.abs(ax - 0.0)))
    le_x = float(pts[np.argmin(np.abs(x_proj - (x_proj.min() if le_at_low else x_proj.max())))][0])

    le_i, _ = ForceProjector._section_ends(pts, c_hat, SPAN)
    ends_le_at_hi = le_i == int(np.argmax(x_proj))
    # both statements in the same frame before comparing them
    deck_le_at_hi = not le_at_low
    agree = deck_le_at_hi == ends_le_at_hi

    print(
        f"{k:3d} {name:>16} {score:8.4f} {le_to_te:+9.1f} {mp_factor:+10.1f} "
        f"{le_x:+8.3f} {'hi-end' if ends_le_at_hi else 'lo-end':>14} {str(agree):>6}"
    )

# estimator convention: a rigid rotation of the tip ring about +z
phys = t._physical_stations(coords)
tip = np.where(np.abs(coords[:, 2] - phys[-1]) < t.STATION_GAP_TOLERANCE)[0]
angle = 0.01  # rad, positive about +z
u = np.zeros(6 * len(coords))
centre = coords[tip].mean(axis=0)
for i in tip:
    u[6 * i : 6 * i + 3] = angle * np.cross(np.array([0.0, 0.0, 1.0]), coords[i] - centre)
omega = float(t._ring_kinematics(coords, u, tip)["omega"])
print(
    f"\nestimator: a rigid rotation of the tip ring by +{angle} rad about +z "
    f"reads omega = {omega:+.6f} rad -> "
    f"{'positive' if omega > 0 else 'negative'}; the module calls omega < 0 nose-down, "
    f"and with the leading edge at +x a +z rotation moves it toward +y (downwind)"
)
