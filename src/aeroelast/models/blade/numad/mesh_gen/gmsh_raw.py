"""Raw-mesh container and chart helpers shared by the gmsh skin mesher and its material mapping.

Moved out of ``windio2_mesh.py`` (odd/tasks/gmsh-official-mesher.md T10, S9) so the gmsh mesher no longer depends on the
windio2 mesher; ``windio2_mesh`` re-exports ``RawMesh``, ``_centroid`` and ``Pchip`` until it is deleted.
"""
from dataclasses import dataclass

import numpy as np

from aeroelast.models.blade.numad.objects.windio2_precision import Pchip  # noqa: F401  (re-exported)


def _centroid(parameters, elements):
    """Mean chart coordinate over each element's true corners (triangles repeat one)."""
    live = (elements >= 0)[:, :, None]
    return (parameters*live).sum(axis=1)/live.sum(axis=1)


@dataclass(frozen=True)
class RawMesh:
    nodes: np.ndarray
    elements: np.ndarray
    cells: tuple
    surfaces: tuple
    element_parameters: np.ndarray
    bases: np.ndarray  # rows: geometric shell x, y, normal; NOT source material axes
    attachments: dict
    material_stacks: tuple | None = None

    @property
    def material_angles(self):
        """Element angles in mapped bottom-to-top order, or None if not requested."""
        if self.material_stacks is None:
            return None
        return tuple(tuple(layer.angle_deg for layer in stack.layers)
                     for stack in self.material_stacks)
