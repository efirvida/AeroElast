"""D-Tube cantilever (Camarena & Anderson 2025, Section 3.1.1) against CalculiX S8R.

The live test (``test_dtube_camarena_bending.py``) finds the MITC4 tip deflection 1.07 % below the
prismatic-beam value at its finest mesh. This file puts a second shell code on the same tube
(semicircle of mid-surface radius R below the web line, flat web on the diameter, wall t = 80 mm,
steel, clamped at z = 0, 63 kPa inward normal pressure on the semicircle) so the residual can be
attributed: either both shells sit away from the beam value (a beam-theory effect), or only ours
does (an element effect).
"""

from __future__ import annotations

import math

import pytest

from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_relative_error
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck

R, T, L = 2.25, 0.08, 100.0
E, NU, P = 2.0e11, 0.3, 63.0e3

A_SEC = math.pi * R * T + 2.0 * R * T
Y_C = -(2.0 * R / math.pi) * (math.pi * R * T) / A_SEC
I_C = 0.5 * math.pi * R**3 * T - A_SEC * Y_C**2
DELTA_BEAM = P * 2.0 * R * L**4 / (8.0 * E * I_C)


def _dtube_deflection(tmp_path, ccx_bin: str, n_arc: int, n_web: int, n_z: int) -> float:
    n_u = n_arc + n_web

    def point(u: float, v: float):
        s = u * n_u
        z = v * L
        if s <= n_arc:
            th = math.pi * s / n_arc
            return (R * math.cos(th), -R * math.sin(th), z)
        return (-R + 2.0 * R * (s - n_arc) / n_web, 0.0, z)

    inp = tmp_path / f"dtube_{n_arc}_{n_z}.inp"
    coords = write_s8r_surface_deck(
        inp,
        point,
        n_u,
        n_z,
        young=E,
        poisson=NU,
        thickness=T,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(z) < 1e-9 else [],
        periodic_u=True,
        pressure=(lambda ei, ej: ei < n_arc, P),
        print_nodes=[
            node_id(i, 2 * n_z, n_u, periodic=True) for i in range(0, 2 * n_arc + 1)
        ],
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    arc_tip = [node_id(i, 2 * n_z, n_u, periodic=True) for i in range(0, 2 * n_arc + 1)]
    del coords
    return abs(sum(float(disp[n][1]) for n in arc_tip) / len(arc_tip))


@pytest.mark.parametrize("n_arc, n_web, n_z", [(32, 16, 100)])
def test_dtube_tip_deflection_calculix_s8r_against_beam_theory(tmp_path, n_arc, n_web, n_z):
    tip = _dtube_deflection(tmp_path, ccx_bin_or_skip(), n_arc, n_web, n_z)
    assert_relative_error(
        tip,
        DELTA_BEAM,
        tol=0.05,
        reference_name="prismatic cantilever w L^4 / (8 E I_c) about the section centroid",
        kind="analytical",
        what=f"D-Tube tip deflection, CalculiX S8R {n_arc}+{n_web} x {n_z}",
    )
