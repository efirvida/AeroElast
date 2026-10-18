"""``classify_modes`` (tools/run_s4_rotating_modal.py) judges a mode by the whole blade, not one node.

T12 of odd/tasks/beamdyn-arbiter-midsurface.md: at es 0.25 (oml, 7.56 rpm) a mode with a small
tip amplitude was labelled "2nd flap" because the label came from the single highest-z node.
These cases use a synthetic clamped thin-walled tube (chord along x, thickness along y, span
along +z, the blade convention) with shapes whose family is known by construction
(odd/tasks/s4-mode-classifier.md, T2).
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

from run_s4_rotating_modal import classify_modes  # noqa: E402

DPN = 6
SPAN, N_RINGS, N_CIRC = 100.0, 121, 32


def _tube():
    z = np.repeat(np.linspace(0.0, SPAN, N_RINGS), N_CIRC)
    th = np.tile(np.linspace(0.0, 2.0 * np.pi, N_CIRC, endpoint=False), N_RINGS)
    xyz = np.c_[2.0 * np.cos(th), 0.5 * np.sin(th), z]
    nodes = [SimpleNamespace(id=i, coords=tuple(p)) for i, p in enumerate(xyz)]
    mesh = SimpleNamespace(nodes=nodes, node_id_to_index={i: i for i in range(len(nodes))})
    free = np.flatnonzero(np.repeat(z > 0.0, DPN))  # root ring clamped
    return mesh, xyz, th, free


def _modes(shapes, free):
    """Stack (n_nodes, 3) translation fields into the solver's (n_modes, n_free) layout."""
    out = np.zeros((len(shapes), free.size))
    for k, u in enumerate(shapes):
        full = np.zeros(u.shape[0] * DPN)
        for d in range(3):
            full[d::DPN] = u[:, d]
        out[k] = full[free]
    return out


def _bend(xyz, ux, uy):
    s = xyz[:, 2] / SPAN
    w = s**2 * (3.0 - s)  # cantilever-like profile
    return np.c_[ux * w, uy * w, np.zeros_like(s)]


def test_local_wall_mode_whose_tip_node_moves_in_y_is_not_flap():
    mesh, xyz, th, free = _tube()
    s = xyz[:, 2] / SPAN
    u = np.zeros_like(xyz)
    # Wall ripple near the tip: cos(4 theta) has no net section translation or rotation,
    # but at the highest-z node chosen by argmax (theta = 0) it moves purely in y.
    u[:, 1] = np.where(s > 0.85, np.cos(4.0 * th) * (s - 0.85), 0.0)
    assert classify_modes(mesh, _modes([u], free), free) == ["local"]


def test_pure_torsion_is_torsion_not_flap():
    mesh, xyz, _, free = _tube()
    s = xyz[:, 2] / SPAN
    u = np.c_[-xyz[:, 1] * s, xyz[:, 0] * s, np.zeros_like(s)]
    assert classify_modes(mesh, _modes([u], free), free) == ["torsion"]


def test_coupled_bending_keeps_its_family():
    """Previous behaviour still holds: IEA-15 1F has tip |ux|/|uy| ~ 0.25, 1E |uy|/|ux| ~ 0.34."""
    mesh, xyz, _, free = _tube()
    shapes = [_bend(xyz, 0.25, 1.0), _bend(xyz, 1.0, 0.34), _bend(xyz, 1.0, 1.0)]
    assert classify_modes(mesh, _modes(shapes, free), free) == ["flap", "edge", "mixed"]
