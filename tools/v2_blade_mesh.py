"""The one blade-mesh constructor for the tools: gmsh skin mesher on the official windIO 2.x inputs.

gmsh-official-mesher T15 (S9): the NuMAD ``BladeMesh`` on windIO 1.x yaml / UTD xlsx inputs is retired; every tool that
needs a blade shell mesh builds it here.  Default sizing is the reference mesh of AGENTS.md rule 9 (chord-scaled N20,
span aspect 1); ``element_size`` selects the uniform-size path and must stay in the certified range [0.05, 0.25] m
(``element_size >= 0.3`` is rejected, ``docs/blade_mesh_size_rule.md``).
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
V2_IEA15 = REPO_ROOT / "tests" / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"
V2_IEA34 = REPO_ROOT / "tests" / "IEA34" / "reference" / "upstream" / "IEA-3.4-130-RWT.yaml"
REFERENCE_CHORD_DIVISIONS = 20
REFERENCE_SPAN_ASPECT = 1.0


def blade_generator(yaml_file=V2_IEA15, element_size: float | None = None, skin_reference: str = "mid"):
    """``GmshBladeMesh`` on a windIO 2.x yaml: reference mesh N20 AR 1 unless ``element_size`` is given."""
    from aeroelast.core.mesh.gmsh_blade import GmshBladeMesh

    if element_size is not None and element_size >= 0.3:
        raise ValueError(f"element_size {element_size} >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")
    sizing = ({"element_size": element_size} if element_size is not None
              else {"chord_divisions": REFERENCE_CHORD_DIVISIONS, "span_aspect": REFERENCE_SPAN_ASPECT})
    return GmshBladeMesh(yaml_file=str(yaml_file), skin_reference=skin_reference, **sizing)
