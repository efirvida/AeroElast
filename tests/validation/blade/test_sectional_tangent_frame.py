"""S-1 bending stiffness read in the section tangent frame (odd/tasks/l3-physical-validation.md L28).

``SectionalExtractor.ei_from_static_response`` reads moment and curvature in global X/Y, while the BeamDyn
6x6 is in the local tangent frame. On the prebent IEA-15 blade a global-Y moment has a torsion component
about the tilted beam axis; the blade is torsionally soft (EI/GJ ~ 46-69), so
1/EI_global = cos^2(phi)/EI_local + sin^2(phi)/GJ and EI_edge reads low (L27: 3.6-6 pp).

Discriminant (cheapest prebent prismatic segment of the L27 study, ``logs/check/l3-prismatic``): a 24 m
segment with every geometric and layup grid held at r = 62.91 m, es 0.25, skin "mid", reference axis
straight (A) or with the real prebend slope (F1, 0.036 rad, ~2 deg; the same tilt as the real blade
mid-span). The reference is the PURE-MOMENT edge/flap EI of the STRAIGHT segment A (tip force couple,
rigid-section rotation slope), a quantity the prebend cannot change: the section is identical.

Marks: slow (two 24 m meshes, ~6k nodes each, a few solves; about a minute).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")

REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_YAML = REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"
KEYS_MOD = {"offset_y_pa", "thickness", "fiber_orientation", "width", "start_nd_arc", "end_nd_arc"}
S_MID, S_A, S_B, LSEG, NODES, ES = 62.91 / 117.0, 51.38 / 117.0, 74.67 / 117.0, 24.0, 25, 0.25
WINDOW = (0.3, 0.7)  # fraction of the segment length; clear of the clamp and the tip load introduction

# Tolerances, from the measured L28 scatter (SLURM 11622461, logs/check/l3-s1-frame/measnew-11622461.log).
# Tangent-frame EI_edge over the 0.3-0.7 L stations, relative to the straight pure-moment EI: straight
# -0.78 % (stations -1.5..-0.1 %), prebent -1.28 % (stations -1.9..-1.0 %); EI_flap -1.04 % / -1.06 %.
# The ~1 % offset is tip-force vs pure-moment loading. TOL_TANGENT is 1.3x the worst single station;
# TOL_PREBEND bounds the straight-vs-prebent difference of the tangent reading (measured 0.5 pp).
# The default global reading misses by -12.5 % (EI_edge), 5x TOL_TANGENT (asserted: beyond 4x).
TOL_TANGENT = 0.025
TOL_PREBEND = 0.015


def _at(gv, s, full, modified):
    from scipy.interpolate import PchipInterpolator

    g, v = np.asarray(gv["grid"], float), np.asarray(gv["values"], float)
    if not modified:
        return float(np.interp(s, g, v))
    vals = np.zeros(len(full))
    m = (full >= g[0]) & (full <= g[-1])
    vals[m] = PchipInterpolator(g, v, extrapolate=True)(full[m])
    return float(np.interp(s, full, vals))


def write_prismatic_yaml(out: Path, prebend: float) -> None:
    """Official yaml cropped to a 24 m prismatic segment (r = 62.91 m values), prebend slope x ``prebend``.

    Same construction as logs/check/l3-prismatic/flatten_crop.py (A: ``prebend=0``, F1: ``prebend=1``).
    """
    import yaml

    d = yaml.load(open(SRC_YAML), Loader=yaml.Loader)
    b = d["components"]["blade"]
    o, st = b["outer_shape_bem"], b["internal_structure_2d_fem"]
    full = np.asarray(o["reference_axis"]["z"]["grid"], float)
    u = np.linspace(0.0, 1.0, NODES)
    s = S_A + u * (S_B - S_A)

    def const(gv, mod):
        return {"grid": [float(x) for x in u], "values": [_at(gv, S_MID, full, mod)] * NODES}

    for k in ("chord", "twist", "pitch_axis"):
        o[k] = const(o[k], False)
    label = o["airfoil_position"]["labels"][int(np.argmin(np.abs(np.asarray(o["airfoil_position"]["grid"]) - S_MID)))]
    o["airfoil_position"] = {"grid": [0.0, 1.0], "labels": [label] * 2}
    xr = o["reference_axis"]["x"]
    x_mid = np.interp(S_MID, xr["grid"], xr["values"])
    ra = {
        "x": {"grid": [float(v) for v in u],
              "values": [float(np.interp(x, xr["grid"], xr["values"]) - x_mid) * prebend for x in s]},
        "y": {"grid": [0.0, 1.0], "values": [0.0, 0.0]},
        "z": {"grid": [float(v) for v in u], "values": [LSEG * float(v) for v in u]},
    }
    o["reference_axis"] = ra
    st["reference_axis"] = ra
    for part in ("layers", "webs"):
        for layer in st[part]:
            for key, gv in list(layer.items()):
                if isinstance(gv, dict) and "grid" in gv and "values" in gv:
                    layer[key] = const(gv, key in KEYS_MOD)
                elif isinstance(gv, dict) and "fixed" in gv:
                    layer[key] = {"fixed": gv["fixed"], **const({"grid": gv["grid"], "values": gv["values"]}, False)}
    yaml.dump(d, open(out, "w"), Dumper=yaml.Dumper, default_flow_style=None, sort_keys=False)


def build(yaml_file: Path):
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.postprocess.sectional import SectionalExtractor

    gen = BladeMesh(yaml_file=str(yaml_file), element_size=ES, skin_reference="mid")
    mesh = gen.generate(renumber="rcm", verbose=False)
    props = build_rust_properties(gen.numad_mesh_data)
    return mesh, props, SectionalExtractor(mesh, props)


def pure_moment_ei(mesh, ext):
    """Pure tip-couple EI(z) (flap, edge) of a segment: ring-fit section rotation slope, 1/eig(compliance).

    Zero-net-force nodal couple (Mx, My) on the tip ring, root clamped (6 dof); the section rotation of
    every 0.2 m node ring is the plane fit of u_z, its slope over +-2.5 m is the curvature.
    """
    from scipy.sparse.linalg import splu

    coords = np.asarray(mesh.coords_array, float)
    nn, n_dof = len(coords), ext.K.shape[0]
    root = np.array(sorted(mesh.node_id_to_index[n] for n in mesh.get_node_set("RootNodes").nodes))
    order = np.argsort(coords[:, 2])
    cut = np.where(np.diff(coords[order, 2]) > 0.08)[0] + 1
    rings = [order[a:b] for a, b in zip(np.r_[0, cut], np.r_[cut, nn], strict=True)]
    rz = np.array([coords[g, 2].mean() for g in rings])
    tip = rings[-1]
    free = np.setdiff1d(np.arange(n_dof), (root[:, None] * 6 + np.arange(6)).ravel())
    lu = splu(ext.K[free][:, free].tocsc())
    xc, yc = (coords[tip, i] - coords[tip, i].mean() for i in (0, 1))
    gram = np.array([[(xc * xc).sum(), (xc * yc).sum()], [(xc * yc).sum(), (yc * yc).sum()]])
    m_unit = 1e6
    comp = np.zeros((len(rz), 2, 2))
    for j, (mx, my) in enumerate(((m_unit, 0.0), (0.0, m_unit))):
        a, b = np.linalg.solve(gram, np.array([-my, mx]))
        f = np.zeros(n_dof)
        f[tip * 6 + 2] = a * xc + b * yc
        u = np.zeros(n_dof)
        u[free] = lu.solve(f[free])
        disp = u.reshape(nn, 6)[:, :3]
        th = []
        for g in rings:
            basis = np.stack([np.ones(len(g)), coords[g, 0], coords[g, 1]], 1)
            c = np.linalg.lstsq(basis, disp[g, 2], rcond=None)[0]
            th.append([c[2], -c[1]])
        th = np.array(th)
        for i, z in enumerate(rz):
            m = np.abs(rz - z) <= 2.5
            for k in range(2):
                comp[i, k, j] = np.polyfit(rz[m], th[m, k], 1)[0] / m_unit
    ev = np.array([np.sort(1.0 / np.linalg.eigvalsh(0.5 * (c + c.T))) for c in comp])
    return rz, ev[:, 0], ev[:, 1]


def window_mean(z, values, length):
    m = (z > WINDOW[0] * length) & (z < WINDOW[1] * length)
    return float(np.nanmean(values[m]))


@pytest.fixture(scope="module")
def readings(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("tangent_frame")
    out = {}
    for tag, slope in (("straight", 0.0), ("prebent", 1.0)):
        yml = tmp / f"{tag}.yaml"
        write_prismatic_yaml(yml, slope)
        mesh, props, ext = build(yml)
        length = float(np.asarray(mesh.coords_array)[:, 2].max())
        res = {"L": length}
        res["global"] = ext.ei_from_static_response(props)
        res["tangent"] = ext.ei_from_static_response(props, frame="tangent")
        if tag == "straight":
            res["pure"] = pure_moment_ei(mesh, ext)
            objects = (mesh, props, ext)
        out[tag] = res
    out["objects"] = objects
    return out


@pytest.fixture(scope="module")
def readings_objects(readings):
    return readings["objects"]


def legacy_global_ei(mesh, props, ext, load=1e6):
    """Independent copy of the pre-L28 global-frame algorithm (deflection spline), to pin the default path."""
    from scipy.interpolate import UnivariateSpline

    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    model = {"solver": {}, "elements": {"element_family": ElementFamily.SHELL,
                                        "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    membership = ext.station_membership()
    sts = sorted(membership)
    zs = np.array([membership[st]["z"] for st in sts])
    tip_nodes = membership[sts[-1]]["nodes"]
    idx = mesh.node_id_to_index
    rootd = sorted(idx[n] * 6 + d for n in mesh.get_node_set("RootNodes").nodes for d in range(6))
    out = {}
    for name, comp in (("flap", 1), ("edge", 0)):
        s = StaticLinearSolver(mesh, model)
        s.add_dirichlet_conditions([DirichletCondition(rootd, 0.0)])
        s.add_nodal_loads([NodalLoad([int(t) * 6 + comp], [load / len(tip_nodes)]) for t in tip_nodes])
        u = s.solve()
        prof = np.array([u[membership[st]["nodes"] * 6 + comp].mean() for st in sts])
        kappa = np.abs(UnivariateSpline(zs, prof, s=1e-2).derivative(2)(zs))
        out[name] = np.where(kappa > 1e-12, load * (zs.max() - zs) / np.maximum(kappa, 1e-12), np.nan)
    return zs, out["flap"], out["edge"]


@pytest.mark.slow
class TestTangentFrameExtraction:
    def test_default_frame_misses_the_straight_segment_ei(self, readings):
        """RED before the fix of L28: the default (global) reading of the prebent segment is far below the
        straight pure-moment EI_edge (-12.5 % measured; -10.1 % on the straight segment itself, the short-segment
        smoothing bias of the deflection-spline estimator, plus -2.4 pp from the prebend torsion coupling).
        """
        length = readings["straight"]["L"]
        pure = readings["straight"]["pure"]
        ref = window_mean(pure[0], pure[2], length)
        g_p = window_mean(readings["prebent"]["global"][0], readings["prebent"]["global"][2], length)
        assert g_p / ref - 1 < -4 * TOL_TANGENT, f"global/pure = {g_p / ref - 1:+.2%}"

    def test_tangent_frame_recovers_straight_segment_edge_ei(self, readings):
        length = readings["straight"]["L"]
        pure = readings["straight"]["pure"]
        ref = window_mean(pure[0], pure[2], length)
        z, _, edge = readings["prebent"]["tangent"]
        got = window_mean(z, edge, length)
        assert abs(got / ref - 1) < TOL_TANGENT, f"EI_edge tangent/pure-straight = {got / ref - 1:+.2%}"

    def test_tangent_frame_recovers_straight_segment_flap_ei(self, readings):
        length = readings["straight"]["L"]
        pure = readings["straight"]["pure"]
        ref = window_mean(pure[0], pure[1], length)
        z, flap, _ = readings["prebent"]["tangent"]
        got = window_mean(z, flap, length)
        assert abs(got / ref - 1) < TOL_TANGENT, f"EI_flap tangent/pure-straight = {got / ref - 1:+.2%}"

    def test_tangent_frame_is_insensitive_to_the_prebend(self, readings):
        """Same section, tilted axis: the tangent-frame reading must not move (measured 0.5 pp)."""
        length = readings["straight"]["L"]
        for k, name in ((1, "flap"), (2, "edge")):
            s = window_mean(readings["straight"]["tangent"][0], readings["straight"]["tangent"][k], length)
            p = window_mean(readings["prebent"]["tangent"][0], readings["prebent"]["tangent"][k], length)
            assert abs(p / s - 1) < TOL_PREBEND, f"{name}: prebent/straight = {p / s - 1:+.2%}"

    def test_default_path_is_unchanged(self, readings_objects):
        """frame="global" (the default) reproduces the pre-L28 algorithm bit for bit on the straight segment."""
        mesh, props, ext = readings_objects
        new = ext.ei_from_static_response(props)
        old = legacy_global_ei(mesh, props, ext)
        for a, b in zip(new, old, strict=True):
            np.testing.assert_array_equal(a, b)
        explicit = ext.ei_from_static_response(props, frame="global")
        for a, b in zip(new, explicit, strict=True):
            np.testing.assert_array_equal(a, b)

    def test_unknown_frame_is_rejected(self, readings_objects):
        _, props, ext = readings_objects
        with pytest.raises(ValueError, match="frame must be"):
            ext.ei_from_static_response(props, frame="local")
