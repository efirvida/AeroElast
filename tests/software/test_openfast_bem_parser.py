"""The AeroDyn driver reader returns the last-revolution mean of every channel, not one azimuth.

Thrust and torque were already revolution means; the per-node angle of attack and force coefficients were
read from the last row, i.e. one azimuth chosen by where ``tmax`` fell. With the 6 deg shaft tilt of the IEA
15 MW rotor alpha varies 2-7 deg over a revolution (measured, logs/check/l3/close/l4_azimuth_check.log), so
the L4 spanwise rows moved 5-7 pp of max|alpha| when dt changed. CCBlade integrates over azimuth, so the
reference must be the azimuth mean.
"""

import numpy as np

from tests.support.openfast_bem import parse_aerodyn_out

NODES = 3
CHANNELS = ("Alpha", "Cl", "Cd", "Cn", "Ct", "Fn", "Ft", "AxInd", "TnInd")


def _write_out(path, revolutions=3, step_deg=10.0):
    header = ["Time", "Case", "Azimuth", "HWindSpeedX", "RotSpeed", "BldPitch1",
              "RtFldFxh", "RtFldMxh", "RtFldCp", "RtFldCt"]
    header += [f"AB1N{i:03d}{c}" for i in range(1, NODES + 1) for c in CHANNELS]
    az = np.tile(np.arange(0.0, 360.0, step_deg), revolutions)
    rows = []
    for k, a in enumerate(az):
        s = np.sin(np.deg2rad(a))
        row = [float(k), 1.0, a, 9.0, 7.56, 0.0, 1000.0 + 10.0 * s, 500.0 + 5.0 * s, 0.4, 0.8]
        for i in range(NODES):
            for c in CHANNELS:
                row.append(10.0 * (i + 1) + 2.0 * s + (0.5 if c == "Cn" else 0.0))
        rows.append(row)
    lines = [" " + "  ".join(header), " " + "  ".join("(-)" for _ in header)]
    lines += [" " + "  ".join(f"{v:.8e}" for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return az


def test_per_node_channels_are_the_last_revolution_mean_not_the_last_row(tmp_path):
    out = tmp_path / "case.out"
    az = _write_out(out)
    assert np.sin(np.deg2rad(az[-1])) != 0.0  # the last row is off the mean: the test can fail
    ref = parse_aerodyn_out(out)[0]
    # alpha = 10 (i+1) + 2 sin(az); the mean of sin over a full revolution is 0
    np.testing.assert_allclose(ref.alpha_deg, [10.0, 20.0, 30.0], atol=1e-6)
    np.testing.assert_allclose(ref.cn, [10.5, 20.5, 30.5], atol=1e-6)
    np.testing.assert_allclose(ref.fn, [10.0, 20.0, 30.0], atol=1e-6)


def test_thrust_and_torque_keep_their_revolution_mean(tmp_path):
    out = tmp_path / "case.out"
    _write_out(out)
    ref = parse_aerodyn_out(out)[0]
    assert abs(ref.thrust - 1000.0) < 1e-6 and abs(ref.torque - 500.0) < 1e-6


def test_time_azimuth_and_operating_point_still_come_from_the_last_row(tmp_path):
    out = tmp_path / "case.out"
    _write_out(out)
    ref = parse_aerodyn_out(out)[0]
    assert ref.time == 3 * 36 - 1 and ref.wind_speed == 9.0 and ref.rpm == 7.56
