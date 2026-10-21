"""Frequency extraction of `tools/bd_rotating_freevib.py` on synthetic BeamDyn-like signals.

The driver's tip response is a constant, a forced 1P term (the global edge load turning in the
rotating frame) and free modes; the fit must return the modes, never the 1P line.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bd_rotating_freevib import fit_frequencies  # noqa: E402


def _signals(t, one_p_hz, modes):
    flap = 0.02 + 0.3 * np.cos(2 * np.pi * one_p_hz * t) - 0.1 * np.sin(2 * np.pi * one_p_hz * t)
    edge = -0.01 + 0.5 * np.sin(2 * np.pi * one_p_hz * t)
    for f, sigma, a_flap, a_edge in modes:
        decay = np.exp(-sigma * t)
        flap = flap + a_flap * decay * np.cos(2 * np.pi * f * t + 0.3)
        edge = edge + a_edge * decay * np.sin(2 * np.pi * f * t - 0.7)
    return np.vstack([flap, edge])


def test_fit_recovers_modes_and_ignores_the_1p_line():
    t = np.arange(0.0, 200.0, 0.005)
    modes = [(0.5603, 1e-3, 1.0, 0.05), (0.8421, 5e-4, 0.03, 0.8), (1.507, 2e-3, 0.2, 0.01)]
    fitted = fit_frequencies(t, _signals(t, 1 / 6, modes), one_p_hz=1 / 6, band=(0.2, 3.0))
    got = sorted(m["f_hz"] for m in fitted["modes"])
    assert len(got) == 3
    np.testing.assert_allclose(got, [m[0] for m in modes], rtol=1e-7)
    assert fitted["relative_residual"] < 1e-6


def test_fit_never_leaves_the_band_for_a_strong_out_of_band_mode():
    # Box parked (SLURM 11621046): an unbounded fit wandered to f = -20.4 Hz (a torsion-range line).
    t = np.arange(0.0, 100.0, 0.002)
    modes = [(0.8218, 0.0, 0.01, 1.0), (1.2067, 0.0, 1.0, 0.01), (20.4, 1e-3, 3.0, 3.0)]
    fitted = fit_frequencies(t, _signals(t, 0.0, modes), one_p_hz=0.0, band=(0.2, 6.0))
    got = [m["f_hz"] for m in fitted["modes"]]
    assert all(0.2 < f < 6.0 for f in got)
    np.testing.assert_allclose(sorted(got), [0.8218, 1.2067], rtol=1e-6)


def test_fit_reports_the_dominant_channel_of_each_mode():
    t = np.arange(0.0, 200.0, 0.005)
    modes = [(0.5603, 0.0, 1.0, 0.05), (0.8421, 0.0, 0.03, 0.8)]
    fitted = fit_frequencies(t, _signals(t, 1 / 6, modes), one_p_hz=1 / 6, band=(0.2, 3.0))
    by_f = {round(m["f_hz"], 4): m["amplitudes"] for m in fitted["modes"]}
    assert by_f[0.5603][0] > by_f[0.5603][1]
    assert by_f[0.8421][1] > by_f[0.8421][0]
