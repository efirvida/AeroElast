"""IEA-22-280-RWT (windIO 2.0): the v2 reader must not depend on IEA-15 layer names."""
import hashlib
import json
from pathlib import Path

import pytest

from aeroelast.models.blade.numad.objects.blade import Blade

ROOT = Path(__file__).resolve().parents[2] / "IEA22MW" / "reference" / "upstream"
IEA22 = ROOT / "IEA-22-280-RWT.yaml"


def test_iea22_input_is_the_pinned_upstream_file():
    pin = json.loads((ROOT / "IEA-22-280-RWT.source.json").read_text())
    data = IEA22.read_bytes()
    assert hashlib.sha256(data).hexdigest() == pin["sha256"]
    assert len(data) == pin["bytes"]
    blob = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
    assert blob == pin["git_blob"]


@pytest.fixture(scope="module")
def iea22():
    before = IEA22.read_bytes()
    blade = Blade().read_yaml(str(IEA22))
    assert IEA22.read_bytes() == before
    return blade


@pytest.mark.parametrize("eta, webs", [(0.1, 3), (0.5, 3), (0.9, 2)])
def test_iea22_topology_follows_the_source_webs(iea22, eta, webs):
    section = iea22.windio2_topology.section(eta)
    assert len(section.webs) == webs
    # every web contributes its two shell attachment points to the skin boundary
    assert len(section.points) >= 2 * webs + 2
    assert all(a.nd_arc < b.nd_arc for a, b in section.intervals)


def test_v1_input_without_a_uv_layer_is_still_rejected(tmp_path):
    """The uv role is optional only for windIO 2.0; v1 behaviour is unchanged."""
    text = (Path(__file__).resolve().parents[2] / "IEA-15-240-RWT.yaml").read_text()
    assert text.count("UV_protection") > 0
    renamed = tmp_path / "no_uv.yaml"
    renamed.write_text(text.replace("UV_protection", "gelcoat"))
    with pytest.raises(ValueError, match="No uv components found"):
        Blade().read_yaml(str(renamed))
