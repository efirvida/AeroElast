"""S-1 counterfactual: do the S-1 sectional reds move when the skin sits on the mid-surface?

E8 (odd/tasks/s7-torsion-ratio.md) found the shell skin nodes on the OML with the laminate
stacked outward, whereas the tools that produce beam 6x6 tables stack plies inward from the OML.
This tool measures, with the S-1 test's own mesh and extractors, what that placement alone moves:

  1. S-1 mid-span (0.2-0.85 R) mean |rel err| vs the BeamDyn deck for EI_flap, EI_edge, mass,
     GJ, EA(CLT), on the mesh as built (must reproduce the S-1 test numbers: same code path);
  2. the same five numbers after every skin node is shifted inward by the mean half-thickness of
     its skin elements (the shift of tools/diagnose_oml_reference_surface.py, TE_FLAT left on the
     OML, the variant that inverts no element).

Nothing in the solver, laminates or loads changes. It is a diagnosis, not a fix: a mesh-geometry
change is a model decision for the maintainer.

    scripts/aeroenv.sh python tools/s1_midsurface_counterfactual.py <element_size>

Recorded in odd/tasks/validation-l3-b3-audit.md.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import ea_gj_from_static_response  # noqa: E402

BLADE_YAML = REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"
BEAMDYN = REPO_ROOT / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_BeamDyn_blade.dat"


def inset_to_midsurface(mesh, props) -> float:
    """Shift every outer-skin node inward by the mean half-thickness of its skin elements."""
    X0 = mesh.coords_array.copy()
    idx = mesh.node_id_to_index
    outer = {idx[n] for n in mesh.get_node_set("allOuterShellNods").nodes}
    th_of = {}
    for name in mesh.element_sets_names:
        lam = props.get(name)
        if lam is None or not hasattr(lam, "total_thickness"):
            continue
        for e in mesh.get_element_set(name).elements:
            th_of[getattr(e, "id", e)] = lam.total_thickness
    nsum = defaultdict(lambda: np.zeros(3))
    tsum = defaultdict(float)
    cnt = defaultdict(int)
    for el in mesh.elements:
        ids = [idx[n] for n in el.node_ids] if hasattr(el, "node_ids") else [idx[n.id] for n in el.nodes]
        if not all(a in outer for a in ids) or el.id not in th_of:
            continue
        P = X0[ids]
        nrm = (
            np.cross(P[2] - P[0], P[3] - P[1]) if len(ids) == 4 else np.cross(P[1] - P[0], P[2] - P[0])
        )
        nrm /= np.linalg.norm(nrm)
        for a in ids:
            nsum[a] += nrm
            tsum[a] += th_of[el.id]
            cnt[a] += 1
    zr = np.round(X0[:, 2], 6)
    cz = defaultdict(list)
    for a in outer:
        cz[zr[a]].append(a)
    cent = {z: X0[ns, :2].mean(0) for z, ns in cz.items()}
    votes = sum(np.sign(np.dot(nsum[a][:2], cent[zr[a]] - X0[a, :2])) for a in nsum)
    sign = 1.0 if votes > 0 else -1.0
    te_flat = set()
    te_ids = set()
    for name in mesh.element_sets_names:
        if name.endswith("_TE_FLAT"):
            te_ids.update(getattr(e, "id", e) for e in mesh.get_element_set(name).elements)
    for el in mesh.elements:
        if el.id in te_ids:
            te_flat.update(idx[n] for n in el.node_ids)
    X1 = X0.copy()
    shift = []
    for a in nsum:
        if a in te_flat:
            continue
        n = sign * nsum[a]
        n[2] = 0.0
        n /= np.linalg.norm(n)
        h = 0.5 * tsum[a] / cnt[a]
        X1[a] += h * n
        shift.append(h)
    mesh.coords_array = X1
    print(f"inset: {len(shift)} nodes, sign {sign:+.0f} (votes {votes:+.0f}), "
          f"mean half-thickness {np.mean(shift) * 1e3:.1f} mm, max {np.max(shift) * 1e3:.1f} mm", flush=True)
    return float(np.mean(shift))


def s1_numbers(mesh, props, numad, bd, label) -> dict:
    from aeroelast.postprocess.sectional import SectionalExtractor

    ext = SectionalExtractor(mesh, props)
    sec = ext.extract()
    z, ei_flap, ei_edge = ext.ei_from_static_response(props)
    z_eg, _ea_resp, gj, ea_clt = ea_gj_from_static_response(mesh, props, numad_data=numad)
    mid = (z > 0.2 * z.max()) & (z < 0.85 * z.max())
    midm = (sec.z > 0.2 * sec.z.max()) & (sec.z < 0.85 * sec.z.max())
    mide = (z_eg > 0.2 * z_eg.max()) & (z_eg < 0.85 * z_eg.max())
    out = {
        "EI_flap": np.nanmean(np.abs(ei_flap[mid] / np.interp(z[mid], bd.z, bd.ei_flap) - 1)),
        "EI_edge": np.nanmean(np.abs(ei_edge[mid] / np.interp(z[mid], bd.z, bd.ei_edge) - 1)),
        "mass": np.nanmean(np.abs(sec.mass_dens[midm] / np.interp(sec.z[midm], bd.z, bd.mass_dens) - 1)),
        "GJ": np.nanmean(np.abs(gj[mide] / np.interp(z_eg[mide], bd.z, bd.gj) - 1)),
        "EA_CLT": np.nanmean(np.abs(ea_clt[mide] / np.interp(z_eg[mide], bd.z, bd.ea) - 1)),
    }
    # signed mean too: the sign says whether the shell is stiffer or softer than the deck
    sgn = {
        "EI_flap": np.nanmean(ei_flap[mid] / np.interp(z[mid], bd.z, bd.ei_flap) - 1),
        "EI_edge": np.nanmean(ei_edge[mid] / np.interp(z[mid], bd.z, bd.ei_edge) - 1),
        "mass": np.nanmean(sec.mass_dens[midm] / np.interp(sec.z[midm], bd.z, bd.mass_dens) - 1),
        "GJ": np.nanmean(gj[mide] / np.interp(z_eg[mide], bd.z, bd.gj) - 1),
        "EA_CLT": np.nanmean(ea_clt[mide] / np.interp(z_eg[mide], bd.z, bd.ea) - 1),
    }
    for k in out:
        print(f"RESULT {label} {k}: mean|rel err| {100 * out[k]:.2f} %   mean signed {100 * sgn[k]:+.2f} %", flush=True)
    return out


def main() -> int:
    es = float(sys.argv[1])
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    bd = load_beamdyn_blade(BEAMDYN)
    gen = BladeMesh(yaml_file=str(BLADE_YAML), element_size=es)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    s1_numbers(mesh, props, gen.numad_mesh_data, bd, f"es={es} OML")
    inset_to_midsurface(mesh, props)
    s1_numbers(mesh, props, gen.numad_mesh_data, bd, f"es={es} midsurface")
    return 0


if __name__ == "__main__":
    sys.exit(main())
