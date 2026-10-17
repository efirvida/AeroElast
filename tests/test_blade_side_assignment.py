"""Spar caps and TE reinforcements must land on their own side of the blade.

The WindIO layers are named per side (``*_SS`` suction = low pressure, ``*_PS``
pressure = high pressure). A name-lookup copy/paste once sent both spar caps to
the pressure side and both TE reinforcements to the suction side, which keeps
the total mass right while moving the section stiffness.
"""

from pathlib import Path

import pytest

from aeroelast.models.blade.numad import Blade

YAML = str(Path(__file__).parent / "IEA-15-240-RWT.yaml")


@pytest.fixture(scope="module")
def blade():
    b = Blade()
    b.read_yaml(YAML)
    return b


@pytest.mark.parametrize(
    "ss, ps",
    [("Spar_Cap_SS", "Spar_Cap_PS"), ("TE_reinforcement_SS", "TE_reinforcement_PS")],
)
def test_component_extents_follow_the_side_in_the_name(blade, ss, ps):
    components = blade.definition.components
    assert not components[ss].hpextents and components[ss].lpextents
    assert not components[ps].lpextents and components[ps].hpextents


@pytest.mark.parametrize(
    "name, hp_stacks",
    [
        ("Spar_Cap_SS", False),
        ("Spar_Cap_PS", True),
        ("TE_reinforcement_SS", False),
        ("TE_reinforcement_PS", True),
    ],
)
def test_plies_are_in_the_stacks_of_their_side(blade, name, hp_stacks):
    station = 25
    stacks = blade.stackdb.stacks[:, station]
    on_hp = {j for j, s in enumerate(stacks) if any(g.component == name for g in s.plygroups)}
    assert on_hp, f"{name} has no plies at station {station}"
    hp_ids = {j for j, s in enumerate(stacks) if "_HP_" in s.name}
    lp_ids = {j for j, s in enumerate(stacks) if "_LP_" in s.name}
    assert on_hp <= (hp_ids if hp_stacks else lp_ids)
