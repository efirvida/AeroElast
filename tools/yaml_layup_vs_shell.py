"""Compare standalone YAML integrals with unchanged S-1 extractors at es=0.25 m.

Run via scripts/aeroenv.sh python -B tools/yaml_layup_vs_shell.py.
Writes JSON evidence only under logs/check/layup/. This deliberately invokes
both SectionalExtractor.extract and ea_gj_from_static_response (including their
static solves), rather than silently replacing the published diagnostics.
The per-material mesh integral additionally isolates mapping/geometry from the
nonreciprocal Q convention in the existing S-1 EA diagnostic.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from yaml_layup_independent import load, section, moduli  # noqa: E402
from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import ea_gj_from_static_response  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def main():
    data = load()
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=0.25)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    print("MESH READY", mesh.node_count, flush=True)
    ext = SectionalExtractor(mesh, props)
    sec = ext.extract()
    print("MASS EXTRACTED", flush=True)
    z, _, _, ea = ea_gj_from_static_response(mesh, props, numad_data=gen.numad_mesh_data)
    print("EA EXTRACTED", flush=True)
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    for key, zz, values, reference in (("mass", sec.z, sec.mass_dens, bd.mass_dens), ("EA_CLT", z, ea, bd.ea)):
        mid = (zz > 0.2*zz.max()) & (zz < 0.85*zz.max())
        r = values[mid]/np.interp(zz[mid], bd.z, reference)-1
        print(f"S1_NATIVE {key} signed={100*np.nanmean(r):+.6f}% absolute={100*np.nanmean(abs(r)):.6f}% n={mid.sum()}", flush=True)
    # Same element areas and station containment as S-1, but broken down by
    # mapped material and with physical Q as an additional channel.
    membership = ext.station_membership()
    stations = sorted(membership)
    coords = np.asarray([n.coords for n in mesh.nodes])
    idx = mesh.node_id_to_index
    sections = {s['elementSet']:s for s in gen.numad_mesh_data['sections']}
    mats = {m['name']:m for m in gen.numad_mesh_data['materials']}
    elem_set = {}
    for name, eset in mesh.element_sets.items():
        for e in eset.elements:
            elem_set.setdefault(e.id, name)
    contributions = []
    for st in stations:
        zz = membership[st]['z']
        dz = max(zz-membership[st-1]['z'], 1e-9) if st-1 in membership else 1e-9
        ns = set(membership[st]['nodes'].tolist())
        totals = {}
        for e in mesh.elements:
            name = elem_set.get(e.id)
            if name not in sections or any(idx[n] not in ns for n in e.node_ids):
                continue
            p = coords[[idx[n] for n in e.node_ids]]
            area = (0.5*np.linalg.norm(np.cross(p[2]-p[0], p[3]-p[1])) if len(p) == 4 else 0.5*np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0])))
            for mat_name, h, ang in sections[name]['layup']:
                mat = mats[mat_name]
                elastic = mat['elastic']
                raw = elastic['E']
                orth = isinstance(raw, (list, tuple))
                m = dict(orth=orth, E=raw, G=elastic['G'], nu=elastic['nu'])
                values = totals.setdefault(mat_name, dict(mass=0., ea_q=0., ea_ex=0., ea_s1=0., area=0.))
                volume_per_length = float(h)*area/dz
                for k, v in zip(('ea_q', 'ea_ex', 'ea_s1'), moduli(m, np.deg2rad(ang))*volume_per_length):
                    values[k] += float(v)
                values['mass'] += float(mat['density'] if 'density' in mat else mat['rho'])*volume_per_length
                values['area'] += volume_per_length
        contributions.append(dict(z=zz, materials=totals))
    rows = []
    for eta in np.linspace(.2, .85, 12):
        r = section(data, eta)
        r['shell_ea'] = float(np.interp(eta*117, z, ea))
        r['shell_mass'] = float(np.interp(eta*117, sec.z, sec.mass_dens))
        r['mesh_materials'] = {m: {k: float(np.interp(eta*117, z, [c['materials'].get(m, {}).get(k, 0.) for c in contributions]))
                                      for k in ('mass', 'ea_q', 'ea_ex', 'ea_s1', 'area')} for m in mats}
        for key in ('mass', 'ea_q', 'ea_ex', 'ea_s1'):
            r['mapped_'+key] = sum(v[key] for v in r['mesh_materials'].values())
        rows.append(r)
    print('MATERIAL MAPPING: mean(mesh integral - independent integral), kg/m and GN', flush=True)
    for mat in mats:
        differences = []
        for r in rows:
            independent = {k:sum(l[k] for l in r['layers'] if l['material'] == mat)
                           for k in ('mass', 'ea_q')}
            differences.append([r['mesh_materials'][mat][k]-independent[k] for k in ('mass', 'ea_q')])
        delta = np.mean(differences, axis=0)
        print(f'MATERIAL {mat} mass_delta={delta[0]:+.6f} EA_Q_delta_GN={delta[1]/1e9:+.6f}', flush=True)
    for r in rows:
        print(f"STATION eta={r['eta']:.6f} mass yaml/beam/shell={r['mass']:.5f}/{r['beam_mass']:.5f}/{r['shell_mass']:.5f} "
              f"EA_Q/S1/Ex/beam/shell[GN]={r['ea_q']/1e9:.6f}/{r['ea_s1']/1e9:.6f}/{r['ea_ex']/1e9:.6f}/{r['beam_ea']/1e9:.6f}/{r['shell_ea']/1e9:.6f}", flush=True)
    for label, a, b in (("yaml_mass/beam", 'mass', 'beam_mass'), ("shell_mass/beam", 'shell_mass', 'beam_mass'),
                        ("yaml_mass/shell", 'mass', 'shell_mass'), ("yaml_Q/beam", 'ea_q', 'beam_ea'),
                        ("yaml_S1/beam", 'ea_s1', 'beam_ea'), ("yaml_Ex/beam", 'ea_ex', 'beam_ea'),
                        ("shell_EA/beam", 'shell_ea', 'beam_ea'), ("yaml_Q/shell", 'ea_q', 'shell_ea'),
                        ("yaml_S1/shell", 'ea_s1', 'shell_ea'),
                        ("mapped_Q/beam", 'mapped_ea_q', 'beam_ea'),
                        ("yaml_Q/mapped_Q", 'ea_q', 'mapped_ea_q'),
                        ("mapped_mass/shell_mass", 'mapped_mass', 'shell_mass'),
                        ("mapped_S1/shell_EA", 'mapped_ea_s1', 'shell_ea')):
        ratios = np.array([r[a]/r[b]-1 for r in rows])
        print(f"COMMON_MEAN {label} signed={100*ratios.mean():+.6f}% absolute={100*np.abs(ratios).mean():.6f}%", flush=True)
    output = dict(rows=rows, shell_z=z.tolist(), shell_ea=ea.tolist(), mass_z=sec.z.tolist(),
                  shell_mass=sec.mass_dens.tolist(), mapped_materials=gen.numad_mesh_data['materials'],
                  mesh_material_integrals=contributions)
    path = ROOT / 'logs/check/layup/comparison.json'
    path.write_text(json.dumps(output, indent=2, default=lambda x: np.asarray(x).tolist())+'\n')
    print('WROTE', path, flush=True)


if __name__ == '__main__':
    main()
