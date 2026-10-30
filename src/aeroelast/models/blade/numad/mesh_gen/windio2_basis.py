"""Shell geometry and the explicit S12 SONATA-derived raw adapter policy.

Not a universal windIO stacking rule. Geometry helpers remain independent of
material interpretation; map_stack requires an explicit physical normal.
"""
from dataclasses import dataclass

import numpy as np

from aeroelast.models.blade.numad.objects.windio2_laminate import LayerSample

from .mesh_tools import get_direction_cosines


def shell_bases(nodes, elements):
    """Right-handed rows (x,y,n); reject nonpositive bilinear Jacobians.

    x is the centre's first natural tangent. n follows connectivity. xDir=x,
    xyDir=y reproduces this frame in the existing adapter. A physical fiber's
    positive adapter angle rotates x toward y about n.
    """
    result = []
    for xyz in np.asarray(nodes)[elements]:
        x = xyz[1]-xyz[0]+xyz[2]-xyz[3]
        y = xyz[3]-xyz[0]+xyz[2]-xyz[1]
        if len(np.unique(xyz, axis=0)) != 4 or not np.linalg.norm(np.cross(x, y)) > 0:
            raise ValueError('Nonpositive shell quadrilateral on source geometry')
        if xyz.dtype == np.dtype(np.longdouble):
            # The legacy helper deliberately substitutes global z below an
            # absolute area threshold. A finite v2 band must keep its normal.
            first = x/np.linalg.norm(x)
            normal = np.cross(x, y)
            normal /= np.linalg.norm(normal)
            basis = np.array([first, np.cross(normal, first), normal])
        else:
            basis = get_direction_cosines(x, y)
        local = xyz-xyz[0]
        for r in (-1., 0., 1.):
            for s in (-1., 0., 1.):
                dr = np.array([-(1-s), 1-s, 1+s, -(1+s)])/4
                ds = np.array([-(1-r), -(1+r), 1+r, 1-r])/4
                if not np.dot(np.cross(dr@local, ds@local), basis[2]) > 0:
                    raise ValueError('Nonpositive shell Jacobian on source geometry')
        result.append(basis)
    return np.asarray(result)


def physical_fiber_angle(vector, basis):
    """Degrees of an EXPLICIT physical vector projected into a shell plane.

    This does not interpret windIO fiber_orientation. Passing a source angle
    without constructing its authoritative physical vector is not supported.
    """
    local = np.asarray(basis)[:2] @ np.asarray(vector)
    if not np.isfinite(local).all() or not np.linalg.norm(local) > 0:
        raise ValueError('Physical fiber must have a finite nonzero shell projection')
    return float(np.rad2deg(np.arctan2(local[1], local[0])))


@dataclass(frozen=True)
class MappedLayer:
    source: LayerSample  # unchanged material, thickness and SOURCE angle
    angle_deg: float  # element x toward element y
    fiber_vector: tuple  # physical NuMAD unit vector
    z_bounds_m: tuple  # increasing element-normal distances from reference


@dataclass(frozen=True)
class MappedStack:
    layers: tuple[MappedLayer, ...]  # increasing element z, NOT source order
    reference: str  # outer_mold_line or web_core_center; nodes never relocated
    physical_normal: tuple  # skin outward, or web toward LE

    @property
    def midplane_offset_m(self):
        """Geometric stack centre relative to the unchanged reference surface."""
        if not self.layers:
            return 0.
        return (self.layers[0].z_bounds_m[0] + self.layers[-1].z_bounds_m[1])/2


def map_stack(cell, beam_axis, positive_normal, basis):
    """Map a source cell under the explicitly approved S12 adapter policy.

    Zero is the unit projection of the local beam axis into the shell plane.
    Positive angle rotates toward ``face_normal cross zero``: skin uses its
    outward normal, web LE/TE skins use opposite outward face normals. These
    directed frames are a SONATA-derived adapter convention, NOT an additional
    windIO schema claim (in particular not a universal web angle convention).

    Skin source order goes OML -> inward, with the OML at z=0. Only the source
    web triple [orthotropic LE skin, isotropic core, orthotropic TE skin] is
    supported; its original web chart is the core centre. Core occupies
    [-t_core/2,t_core/2], and skins grow outward from its two faces. No implicit
    centring of an asymmetric stack and no coordinate displacement occurs.
    The caller's right-handed element basis may have either normal direction.
    """
    scalar = np.longdouble if np.asarray(basis).dtype == np.dtype(np.longdouble) else float
    basis = np.asarray(basis, dtype=scalar)
    n = np.asarray(positive_normal, dtype=scalar)
    beam = np.asarray(beam_axis, dtype=scalar)
    if (basis.shape != (3, 3) or not np.isfinite(basis).all()
            or not np.allclose(basis@basis.T, np.eye(3), atol=1e-12, rtol=0)
            or np.dot(np.cross(basis[0], basis[1]), basis[2]) <= 0 or n.shape != (3,)
            or not np.isfinite(n).all() or np.linalg.norm(n) == 0):
        raise ValueError('Material mapping requires a right-handed basis and physical normal')
    n = n/np.linalg.norm(n)
    direction = float(n@basis[2])
    if not np.isclose(abs(direction), 1., atol=1e-12, rtol=0):
        raise ValueError('Physical normal must be parallel to the element normal')
    direction = 1. if direction > 0 else -1.
    if beam.shape != (3,) or not np.isfinite(beam).all() or np.linalg.norm(beam) == 0:
        raise ValueError('Degenerate beam-axis projection into shell plane')
    beam = beam/np.linalg.norm(beam)
    zero = beam - (beam@n)*n
    if np.linalg.norm(zero) <= 1e-12:
        raise ValueError('Degenerate beam-axis projection into shell plane')
    zero /= np.linalg.norm(zero)
    layers = cell.layers
    if cell.surface == 'skin':
        positions, z = [], 0.
        for layer in layers:
            positions.append((z-layer.thickness_m, z, n))
            z -= layer.thickness_m
        reference = 'outer_mold_line'
    else:
        if (len(layers) != 3 or np.shape(layers[0].material['E']) != (3,)
                or np.ndim(layers[1].material['E']) != 0
                or np.shape(layers[2].material['E']) != (3,)):
            raise ValueError('Web material mapping requires one orthotropic LE skin, isotropic core, orthotropic TE skin')
        le, core, te = (layer.thickness_m for layer in layers)
        positions = [(core/2, core/2+le, n), (-core/2, core/2, n),
                     (-core/2-te, -core/2, -n)]
        reference = 'web_core_center'
    mapped = []
    for source, (lo, hi, face) in zip(layers, positions, strict=True):
        angle = np.deg2rad(source.fiber_angle_deg)
        vector = np.cos(angle)*zero + np.sin(angle)*np.cross(face, zero)
        bounds = tuple(sorted((direction*lo, direction*hi)))
        mapped.append(MappedLayer(source, physical_fiber_angle(vector, basis),
                                  tuple(vector), bounds))
    return MappedStack(tuple(sorted(mapped, key=lambda layer: layer.z_bounds_m[0])),
                       reference, tuple(n))
