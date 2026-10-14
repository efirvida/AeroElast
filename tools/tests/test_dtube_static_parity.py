"""Fixed-load D-Tube comparison instrument, independent of FE solvers."""

import math

from tools.diagnose_dtube_static_parity import add_root_reaction_probe, fixed_arc_forces


def test_fixed_arc_forces_keep_both_components_and_resultant():
    loads = fixed_arc_forces(n_arc=32, n_z=100, subdivisions=1, pressure=63e3)
    assert len(loads) == 33 * 101
    fx = sum(force[0] for _, _, force in loads)
    fy = sum(force[1] for _, _, force in loads)
    assert abs(fx) < 1e-6
    assert abs(fy / (63e3 * 4.5 * 100) - 1) < 0.001
    assert any(force[0] != 0.0 for _, _, force in loads)
    assert all(math.isfinite(v) for _, _, force in loads for v in force)


def test_subdividing_the_s8r_grid_preserves_the_load_resultant():
    coarse = fixed_arc_forces(n_arc=32, n_z=100, subdivisions=1, pressure=63e3)
    fine = fixed_arc_forces(n_arc=32, n_z=100, subdivisions=2, pressure=63e3)
    assert len(fine) == 65 * 201
    fcoarse = sum(force[1] for _, _, force in coarse)
    ffine = sum(force[1] for _, _, force in fine)
    assert abs(fcoarse - ffine) / ffine < 0.001


def test_root_reaction_probe_requests_rf_without_changing_the_load(tmp_path):
    deck = tmp_path / "dtube.inp"
    deck.write_text("*BOUNDARY\n1, 1, 1\n*STEP\n*STATIC\n*CLOAD\n1, 2, 5.0\n*END STEP\n")
    add_root_reaction_probe(deck, n_u=2)
    text = deck.read_text()
    assert "*NSET, NSET=ROOT_RF\n1, 2, 3, 4\n*STEP" in text
    assert "*CLOAD\n1, 2, 5.0\n*NODE PRINT, NSET=ROOT_RF\nRF\n*END STEP" in text


def test_root_reaction_probe_preserves_nlgeom_step(tmp_path):
    deck = tmp_path / "dtube_nlgeom.inp"
    deck.write_text("*STEP, NLGEOM, INC=80\n*STATIC\n*END STEP\n")
    add_root_reaction_probe(deck, n_u=1)
    text = deck.read_text()
    assert "*NSET, NSET=ROOT_RF\n1, 2\n*STEP, NLGEOM, INC=80" in text
    assert "*NODE PRINT, NSET=ROOT_RF\nRF\n*END STEP" in text
