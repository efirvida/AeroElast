"""IEA-3.4-130-RWT (windIO 2.1.1): the blade path of the v2 reader accepts 2.0.x / 2.1.x and nothing else.

windIO 2.0.0 -> 2.1.1 changed no blade structure: ``components.blade`` differs only in the wording of the ``twist`` description,
``airfoils`` and ``materials`` are identical, the other differences are control, assembly and the floating platform
(turbine_schema.yaml diffed at tags 2.0.0 / 2.0.1 / 2.1 / 2.1.1, 2026-11-04).
"""
import hashlib
import json
from pathlib import Path

import pytest

from aeroelast.models.blade.numad.objects.blade import Blade
from aeroelast.models.blade.windio import normalize_windio

ROOT = Path(__file__).resolve().parents[2] / "IEA34" / "reference" / "upstream"
IEA34 = ROOT / "IEA-3.4-130-RWT.yaml"


def test_iea34_input_is_the_pinned_upstream_file():
    pin = json.loads((ROOT / "IEA-3.4-130-RWT.source.json").read_text())
    data = IEA34.read_bytes()
    assert hashlib.sha256(data).hexdigest() == pin["sha256"]
    assert len(data) == pin["bytes"]
    assert hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest() == pin["git_blob"]
    assert pin["windIO_version"] == "2.1.1"


@pytest.fixture(scope="module")
def iea34():
    before = IEA34.read_bytes()
    blade = Blade().read_yaml(str(IEA34))
    assert IEA34.read_bytes() == before
    return blade


@pytest.fixture(scope="module")
def iea34_as_20(tmp_path_factory):
    text = IEA34.read_text()
    assert text.startswith("windIO_version: '2.1.1'")
    path = tmp_path_factory.mktemp("iea34") / "IEA-3.4-as-2.0.yaml"
    path.write_text(text.replace("windIO_version: '2.1.1'", "windIO_version: '2.0'", 1))
    return Blade().read_yaml(str(path))


@pytest.mark.parametrize("eta", [0.1, 0.5, 0.9])
def test_iea34_2_1_1_gives_the_topology_of_the_same_file_labelled_2_0(iea34, iea34_as_20, eta):
    a, b = iea34.windio2_topology.section(eta), iea34_as_20.windio2_topology.section(eta)
    assert len(a.webs) >= 2
    assert repr(a.points) == repr(b.points)
    assert repr(a.intervals) == repr(b.intervals)


@pytest.mark.parametrize("version", ["2.2", "3.0", "2"])
def test_other_windio_versions_are_still_rejected(version):
    with pytest.raises(ValueError, match="Unsupported windIO version"):
        normalize_windio({"windIO_version": version})


def test_v1_data_is_returned_unchanged():
    data = {"windIO_version": "1.0", "components": {}}
    assert normalize_windio(data) is data
