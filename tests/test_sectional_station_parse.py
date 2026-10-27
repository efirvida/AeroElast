"""Station membership reads the whole station field of the NuMAD set name (blade input audit T7c).

Set names are ``{segment:02d}_{station:02d}_{label}``; ``:02d`` grows to three digits past station 99,
which ``BladeMesh.max_bay_length`` reaches. Reading ``name[3:5]`` merged station 147 into 14 and broke
every S-1 extraction on the refined mesh (jobs 11623121, 11623128).
"""

from aeroelast.postprocess.sectional import _station_of


def test_three_digit_station():
    assert _station_of("08_147_LP_SPAR") == 147


def test_two_digit_station_unchanged():
    assert _station_of("01_14_HP_TE_REINF") == 14
    assert _station_of("00_05_SW") == 5


def test_non_detailed_sets_are_skipped():
    for name in ("RootNodes", "HP_SPAR", "1_2_X", "ab_12_X", "08_x7_LP"):
        assert _station_of(name) is None
