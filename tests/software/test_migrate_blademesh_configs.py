"""``tools/migrate_blademesh_configs.py`` (gmsh-official-mesher T6b, S6): one pin per mapping rule on small temporary trees.

Each test writes a fake repository under ``tmp_path`` (a windIO 1.x file, the upstream v2 file, one config), runs the tool's
``main`` with ``--write`` and compares the whole output text, so a comment, a key order or a non-generator key that moved fails.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "tools" / "migrate_blademesh_configs.py"
_spec = importlib.util.spec_from_file_location("migrate_blademesh_configs", SCRIPT)
tool = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = tool  # dataclasses resolve their module through sys.modules
_spec.loader.exec_module(tool)

V2_IEA15 = "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
V2_REL = "../IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"  # from tests/case/


@pytest.fixture
def root(tmp_path):
    (tmp_path / "tests/case").mkdir(parents=True)
    (tmp_path / "tests/IEA-15-240-RWT.yaml").write_text("name: v1 turbine\nassembly: {}\n")
    (tmp_path / "tests/IEA-15-240-RWT_rigid.yaml").write_text("name: rigid variant\n")
    for rel, ver in ((V2_IEA15, "2.0"), ("tests/IEA34/reference/upstream/IEA-3.4-130-RWT.yaml", "2.1.1"),
                     ("tests/IEA22MW/reference/upstream/IEA-22-280-RWT.yaml", "2.0")):
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(f"windIO_version: '{ver}'\n")
    return tmp_path


def config(params: str, gtype: str = "BladeMesh", yaml_file: str = "../IEA-15-240-RWT.yaml", tail: str = "",
           blade_file: str = "../IEA-15-240-RWT.yaml") -> str:
    return (
        f"# header comment\nparticipant: Solid\nblade_file: {blade_file}\nmesh:\n  source: generator\n  generator:\n"
        f"    type: {gtype}  # the generator\n    params:\n      yaml_file: {yaml_file}  # blade\n{params}"
        "    coupling_node_set: allOuterShellNods\n  output_file: m.vtu\n"
        "solver:\n  type: LinearDynamicFSIRotorSolver  # keep\n  cfl: 0.5\n" + tail
    )


def migrate(root: Path, text: str, name: str = "x.yaml", write: bool = True, capsys=None):
    path = root / "tests/case" / name
    path.write_text(text)
    args = [str(path), "--root", str(root)] + (["--write"] if write else [])
    assert tool.main(args) == 0
    return path.read_text()


def expected(params: str, gtype: str = "GmshBladeMesh", yaml_file: str = V2_REL, tail: str = "") -> str:
    return config(params, gtype=gtype, yaml_file=yaml_file, tail=tail, blade_file=V2_REL)


LEGACY_COARSE = "0.5"  # a legacy input value the tool must replace; not a mesh this repo builds
SKIN = "      skin_reference: mid  # explicit: what the runner resolved for BladeMesh configs\n"
N20 = "      chord_divisions: 20  # reference mesh N20, span AR 1 (AGENTS.md rule 9)\n      span_aspect: 1\n"


def test_blademesh_type_and_v1_yaml_are_repointed_and_nothing_else_moves(root):
    out = migrate(root, config("      element_size: 0.25\n"))
    assert out == expected(SKIN + "      element_size: 0.25\n")
    assert "blade_file: ../IEA15MW/reference/upstream/IEA-15-240-RWT.yaml\n" in out


def test_rotormesh_keeps_its_type_and_gets_no_skin_reference(root):
    out = migrate(root, config("      element_size: 0.2\n      n_blades: 3\n", gtype="RotorMesh"))
    assert out == expected("      element_size: 0.2\n      n_blades: 3\n", gtype="RotorMesh")


def test_v2_yaml_file_is_left_alone(root):
    text = config("      element_size: 0.25\n", gtype="RotorMesh", yaml_file=V2_REL, blade_file=V2_REL)
    assert migrate(root, text) == text


def test_other_turbines_map_to_their_upstream_v2(root):
    (root / "tests/IEA34").mkdir(exist_ok=True)
    (root / "tests/IEA34/IEA-3.4-130-RWT.yaml").write_text("name: v1 3.4\n")
    out = migrate(root, config("      element_size: 0.25\n", yaml_file="../IEA34/IEA-3.4-130-RWT.yaml"))
    assert out == expected(SKIN + "      element_size: 0.25\n", yaml_file="../IEA34/reference/upstream/IEA-3.4-130-RWT.yaml")


W1 = "tests/archive/windio1_configs"


def test_generator_v1_file_without_v2_twin_is_archived_with_a_readme(root):
    text = config("      element_size: 0.25\n", yaml_file="../IEA-15-240-RWT_rigid.yaml")
    path = root / "tests/case/rigid.yaml"
    path.write_text(text)
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert not path.exists()
    assert (root / W1 / "case/rigid.yaml").read_text() == text
    assert "S8" in (root / W1 / "README.md").read_text()


def test_blade_file_v1_without_v2_twin_is_archived(root):
    text = config("      element_size: 0.25\n", gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file="../IEA-15-240-RWT_rigid.yaml")
    path = root / "tests/case/rigid_bf.yaml"
    path.write_text(text)
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert not path.exists() and (root / W1 / "case/rigid_bf.yaml").read_text() == text


def test_v1_config_without_twin_named_elsewhere_stays_in_place(root, capsys):
    (root / "tools").mkdir()
    (root / "tools/run.sh").write_text("srun fsi tests/case/rigid.yaml\n")
    text = config("      element_size: 0.25\n", yaml_file="../IEA-15-240-RWT_rigid.yaml")
    path = root / "tests/case/rigid.yaml"
    path.write_text(text)
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert path.read_text() == text
    out = capsys.readouterr().out
    assert "SKIPPED" in out and "tools/run.sh" in out


def test_blade_file_alone_is_repointed_on_an_already_migrated_config(root):
    params = "      element_size: 0.25\n"
    out = migrate(root, config(params, gtype="GmshBladeMesh", yaml_file=V2_REL))
    assert out == config(params, gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file=V2_REL)


@pytest.mark.parametrize("quote", ["", '"', "'"])
def test_blade_file_under_aero_is_repointed_keeping_quotes_and_comment(root, quote):
    tail = f"aero:\n  backend: bem\n  blade_file: {quote}../IEA-15-240-RWT.yaml{quote}  # bem input\n  rotor:\n    n_blades: 3\n"
    want = f"aero:\n  backend: bem\n  blade_file: {quote}{V2_REL}{quote}  # bem input\n  rotor:\n    n_blades: 3\n"
    out = migrate(root, config("      element_size: 0.25\n", gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file=V2_REL, tail=tail))
    assert out == config("      element_size: 0.25\n", gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file=V2_REL, tail=want)


def test_blade_file_of_another_turbine_maps_to_its_upstream_v2(root):
    (root / "tests/IEA34").mkdir(exist_ok=True)
    (root / "tests/IEA34/IEA-3.4-130-RWT.yaml").write_text("name: v1 3.4\n")
    out = migrate(root, config("      element_size: 0.25\n", gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file="../IEA34/IEA-3.4-130-RWT.yaml"))
    assert "blade_file: ../IEA34/reference/upstream/IEA-3.4-130-RWT.yaml\n" in out


def test_blade_file_naming_v2_or_a_spreadsheet_is_left_alone(root):
    for bf in (V2_REL, "../blade.xlsx"):
        text = config("      element_size: 0.25\n", gtype="GmshBladeMesh", yaml_file=V2_REL, blade_file=bf)
        assert migrate(root, text) == text


def test_config_without_generator_gets_only_its_blade_file_changed(root):
    text = "participant: Fluid\nblade_file: ../IEA-15-240-RWT.yaml  # c\nbem:\n  yaw: 0.0\n"
    assert migrate(root, text) == f"participant: Fluid\nblade_file: {V2_REL}  # c\nbem:\n  yaw: 0.0\n"


def test_blade_file_migration_is_idempotent(root):
    path = root / "tests/case/x.yaml"
    path.write_text(config("      element_size: 0.25\n", tail="aero:\n  blade_file: ../IEA-15-240-RWT.yaml\n"))
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    once = path.read_text()
    assert once.count(V2_REL) == 3
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert path.read_text() == once
    assert tool.collect([path], root) == []  # nothing left to plan


def test_numad_only_keys_are_dropped(root):
    params = "".join(f"      {k}: 1\n" for k in tool.NUMAD_KEYS) + "      element_size: 0.25\n"
    out = migrate(root, config(params))
    assert out == expected(SKIN + "      element_size: 0.25\n")


@pytest.mark.parametrize("es", ["0.15", "0.2", "0.25"])
def test_element_size_within_0p15_to_0p25_is_kept(root, es):
    out = migrate(root, config(f"      element_size: {es}\n"))
    assert out == expected(SKIN + f"      element_size: {es}\n")


@pytest.mark.parametrize("es", ["0.3", "0.5", "2.0"])
def test_element_size_at_or_above_0p3_becomes_the_n20_reference_in_place(root, es):
    out = migrate(root, config(f"      n_samples: 300\n      element_size: {es}\n      span_grading: chord\n"))
    assert out == expected(SKIN + N20)


def test_absent_element_size_gets_the_n20_reference(root):
    out = migrate(root, config("      n_samples: 300\n"))
    assert out == expected(SKIN + N20)


def test_existing_chord_divisions_are_not_duplicated(root):
    out = migrate(root, config("      chord_divisions: 10\n"))
    assert out == expected(SKIN + "      chord_divisions: 10\n")


def test_element_size_below_0p15_is_kept_and_listed(root, capsys):
    out = migrate(root, config("      element_size: 0.125\n"))
    assert out == expected(SKIN + "      element_size: 0.125\n")
    assert "element_size 0.125 < 0.15 kept" in capsys.readouterr().out


@pytest.mark.parametrize("skin", ["oml", "mid"])
def test_existing_skin_reference_is_never_overridden(root, skin):
    params = f"      skin_reference: {skin}\n      element_size: 0.25\n"
    assert migrate(root, config(params)) == expected(params)


def test_excel_config_is_archived_with_its_subpath_and_a_readme(root):
    path = root / "tests/case/xlsx.yaml"
    path.write_text(config("      excel_file: ../blade.xlsx\n      element_size: 0.25\n"))
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    moved = root / "tests/archive/numad_xlsx_configs/case/xlsx.yaml"
    assert not path.exists()
    assert moved.read_text() == config("      excel_file: ../blade.xlsx\n      element_size: 0.25\n")
    assert "gmsh-official-mesher" in (root / "tests/archive/numad_xlsx_configs/README.md").read_text()


def test_excel_config_named_elsewhere_stays_in_place(root, capsys):
    (root / "tools").mkdir()
    (root / "tools/run.sh").write_text("srun fsi tests/case/xlsx.yaml\n")
    text = config("      excel_file: ../blade.xlsx\n")
    path = root / "tests/case/xlsx.yaml"
    path.write_text(text)
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert path.read_text() == text
    assert "tools/run.sh" in capsys.readouterr().out


def test_dry_run_is_the_default_and_writes_nothing(root, capsys):
    text = config(f"      element_size: {LEGACY_COARSE}\n")
    assert migrate(root, text, write=False) == text
    out = capsys.readouterr().out
    assert "REWRITTEN" in out and "dry-run" in out


def test_second_run_changes_nothing(root):
    path = root / "tests/case/x.yaml"
    path.write_text(config(f"      airfoil_dir: ../a\n      element_size: {LEGACY_COARSE}\n", gtype="RotorMesh"))
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    once = path.read_text()
    assert tool.main([str(path), "--root", str(root), "--write"]) == 0
    assert path.read_text() == once
    plans = tool.collect([path], root)
    assert [p.action for p in plans] == ["unchanged"]


def test_non_generator_keys_are_byte_identical(root):
    text = config(f"      element_size: {LEGACY_COARSE}\n", tail="extra:\n  list: [1, 2]   # flow\n  s: 'q'\n")
    out = migrate(root, text)
    keep = lambda t: [l for l in t.splitlines() if not l.startswith(("      ", "    type", "blade_file"))]  # noqa: E731
    assert keep(out) == keep(text)


@pytest.mark.parametrize("quote", ['"', "'"])
def test_quoted_generator_type_is_migrated_and_keeps_its_quotes(root, quote):
    out = migrate(root, config("      element_size: 0.25\n", gtype=f"{quote}BladeMesh{quote}"))
    assert out == expected(SKIN + "      element_size: 0.25\n", gtype=f"{quote}GmshBladeMesh{quote}")
