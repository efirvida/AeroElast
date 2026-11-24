"""The NuMAD ``--excel-file`` option is gone from the ``aeroelast`` CLI (gmsh-official-mesher T11, S9)."""

import pytest

from aeroelast.cli.aeroelast import main


def test_help_does_not_list_excel_file(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert "--excel-file" not in capsys.readouterr().out


def test_excel_file_option_is_rejected_by_the_parser(capsys, tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["--mesh-generator", "BladeMesh", "--excel-file", str(tmp_path / "blade.xlsx"),
              "--export-mesh", str(tmp_path / "m.vtk")])
    assert exc.value.code == 2
    assert "--excel-file" in capsys.readouterr().err
