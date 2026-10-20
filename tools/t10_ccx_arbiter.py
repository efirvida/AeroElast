#!/usr/bin/env python3
"""T10 same-load shell arbiter (S11, S8, S10); numerical context, NOT a gate.

Predeclared BEFORE runs: |AE/CCX-1| <= 5% for tip-node uy, ring-mean uy,
load work and classified parked frequencies. This is a reporting band, not a
claim of convergence or physical truth: an applicable analytic reference decides.
No paired refinement is performed here. MITC4 versus S8R has documented 4.22%
rated-twist discretisation context at es=0.25; that is not a correction factor.
Load-transfer tolerance: norm(Fdeck-Fae)/norm(Fae) <= 1e-6 (writer rounds to
7 significant digits). Every corner node must match an FRD displacement.
OFFSET probe must expand a 0.01 m laminate on z=0 with +z normal and +0.5
reference offset into z=[-0.01,0] m, at 1e-12 m absolute tolerance.

Usage (SLURM only for blade): scripts/aeroenv.sh python -B
  tools/t10_ccx_arbiter.py {oml,mid,probe} --out logs/check/t10/<run>
The probe tests a one-element composite OFFSET deck. The mid run uses option
(b): independent CCX skin offsets on OML nodes, with centred webs and TE_FLAT.
OFFSET=+0.5 puts the reference nodes on the positive-normal face; thus the
laminate centre moves opposite the element normal. Skin winding is checked across
shared edges; outward sense is determined by signed radial flux, not by assuming
each concave-airfoil normal points away from the node centroid.
Shared skin/web joints and thickness steps differ from averaged-node insetting;
any difference is reported, not treated as exact geometric equivalence.
OpenSeesPy is deliberately excluded: installed build cannot represent +/-45 plies.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'tools'))
CCX = '/scratch/leahk/eduardo.donestevez/venv/bin/ccx'


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def probe(out):
    from tests.support.ccx_io import run_ccx
    deck = out / 'offset_probe.inp'
    deck.write_text('''*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,.5,0,0
6,1,.5,0
7,.5,1,0
8,0,.5,0
*ELEMENT,TYPE=S8R,ELSET=SKIN
1,1,2,3,4,5,6,7,8
*MATERIAL,NAME=MAT
*ELASTIC,TYPE=ENGINEERING CONSTANTS
1e9,1e9,1e9,.3,.3,.3,3.846e8,3.846e8
3.846e8
*DENSITY
1000
*SHELL SECTION,ELSET=SKIN,COMPOSITE,OFFSET=.5
.005,,MAT
.005,,MAT
*BOUNDARY
1,1,6
4,1,6
8,1,6
*STEP
*STATIC
*CLOAD
2,3,1
*NODE FILE,OUTPUT=2D
U
*END STEP
''')
    p = run_ccx(deck, CCX)
    (out / 'offset_probe.stdout').write_text(p.stdout + p.stderr)
    zs = []
    active = False
    for line in deck.with_suffix('.frd').read_text().splitlines():
        if line.startswith('    2C'):
            active = True
        elif active and line.lstrip().startswith('-3'):
            break
        elif active and line.lstrip().startswith('-1'):
            zs.append(float(line[37:49]))
    limits = [min(zs), max(zs)]
    if p.returncode or not np.allclose(limits, [-.01, 0.], rtol=0., atol=1e-12):
        raise RuntimeError(f'OFFSET probe failed: rc={p.returncode}, z={limits}')
    save(out / 'probe.json', {'returncode': p.returncode, 'expanded_z_m': limits,
                              'tail': (p.stdout + p.stderr)[-5000:]})
    print(p.returncode, p.stdout[-4000:])


def metric(xyz, f, u):
    ring = np.flatnonzero(xyz[:, 2] >= xyz[:, 2].max() - 1e-6)
    return {'tip_node_uy_m': float(u[np.argmax(xyz[:, 2]), 1]),
            'tip_ring_uy_m': float(u[ring, 1].mean()),
            'work_J': float(.5 * np.sum(f * u))}


def read_deck_force(deck, n):
    f = np.zeros((n, 6))
    active = False
    for line in deck.read_text().splitlines():
        if line.startswith('**'):
            continue
        if line.startswith('*'):
            active = line.upper().startswith('*CLOAD')
        elif active and line.strip():
            label, dof, value = line.split(',')[:3]
            f[int(label)-1, int(dof)-1] += float(value)
    return f


def frd_modes(path, n):
    """Stream fixed-width DISP blocks; do not use whitespace regex on FRD."""
    blocks = []
    active = False
    with path.open() as handle:
        for line in handle:
            if line.lstrip().startswith('-4') and 'DISP' in line:
                active = True
                u = np.zeros((n, 6))
                seen = set()
            elif active and line.lstrip().startswith('-3'):
                if len(seen) != n:
                    raise RuntimeError(f'modal FRD matched {len(seen)}/{n}')
                blocks.append(u)
                active = False
            elif active and line.lstrip().startswith('-1'):
                k = int(line[3:13])-1
                if 0 <= k < n:
                    u[k, :3] = [float(line[13+12*d:25+12*d]) for d in range(3)]
                    seen.add(k)
    return blocks


def mode_rows(mesh, freqs, modes):
    from diagnose_s4_mode_labels import _label
    tip = np.argmax(mesh.coords_array[:, 2])
    rows = []
    for f, u in zip(freqs, modes):
        p = np.abs(u[:, :3]).sum(axis=0)
        rows.append({'Hz': float(f), 'tip_label': _label(*np.abs(u[tip, :3])),
                     'integrated_label': _label(*p), 'tip_xyz': u[tip, :3].tolist(),
                     'integrated_abs_xyz': p.tolist()})
    return rows


def offset_sections(deck, mesh, props):
    """Post-process only skin section cards; never modify the shared writer."""
    xyz = mesh.coords_array
    idx = mesh.node_id_to_index
    outer = {idx[n] for n in mesh.get_node_set('allOuterShellNods').nodes}
    heights = np.round(xyz[:, 2], 6)
    centers = {z: xyz[[i for i in outer if heights[i] == z], :2].mean(axis=0)
               for z in np.unique(heights[list(outer)])}
    offsets = {}
    flux = 0.
    edges = {}
    skin_names = []
    for name in mesh.element_sets_names:
        if name.endswith('_TE_FLAT') or not hasattr(props.get(name), 'total_thickness'):
            continue
        selected = 0
        for element in mesh.get_element_set(name).elements:
            ids = [idx[n] for n in element.node_ids]
            if not all(i in outer for i in ids):
                continue
            p = xyz[ids]
            normal = (np.cross(p[2]-p[0], p[3]-p[1]) if len(ids) == 4
                      else np.cross(p[1]-p[0], p[2]-p[0]))
            center = np.mean([centers[heights[i]] for i in ids], axis=0)
            flux += float(np.dot(normal[:2], p[:, :2].mean(axis=0)-center))
            selected += 1
            for a, b in zip(ids, ids[1:]+ids[:1]):
                edges.setdefault(tuple(sorted((a, b))), []).append(1 if a < b else -1)
        if selected:
            skin_names.append('E'+name.upper())
    bad = [e for e, signs in edges.items() if len(signs) == 2 and sum(signs) != 0]
    if bad or not flux:
        raise RuntimeError(f'skin winding inconsistent ({len(bad)} edges), flux={flux}')
    offsets = {name: .5*float(np.sign(flux)) for name in skin_names}
    lines = []
    touched = set()
    for line in deck.read_text().splitlines():
        if line.upper().startswith('*SHELL SECTION'):
            fields = dict(token.strip().split('=', 1) for token in line.upper().split(',')
                          if '=' in token)
            region = fields['ELSET']
            if region in offsets:
                line += f', OFFSET={offsets[region]:+.1f}'
                touched.add(region)
        lines.append(line)
    if touched != set(offsets):
        raise RuntimeError('not all skin sections were offset')
    deck.write_text('\n'.join(lines)+'\n')
    return offsets


def blade(ref, out):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('blade runs require SLURM')
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.core.mesh.io.writers import write_ccx_mesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
    from tests.support.ccx_io import run_ccx, parse_frd_disp, parse_ccx_frequencies
    from run_s2_static_cases import run_shell_case

    gen = BladeMesh(yaml_file=str(REPO / 'tests/IEA-15-240-RWT.yaml'),
                    element_size=.25, skin_reference=ref)
    mesh = gen.generate(renumber='rcm')
    props = build_rust_properties(gen.numad_mesh_data)
    xyz = mesh.coords_array
    n = len(xyz)
    ccx_mesh = mesh
    ccx_props = props
    if ref == 'mid':
        oml_gen = BladeMesh(yaml_file=str(REPO / 'tests/IEA-15-240-RWT.yaml'),
                            element_size=.25, skin_reference='oml')
        ccx_mesh = oml_gen.generate(renumber='rcm')
        ccx_props = build_rust_properties(oml_gen.numad_mesh_data)
        conn = lambda m: [[m.node_id_to_index[k] for k in e.node_ids] for e in m.elements]
        if conn(mesh) != conn(ccx_mesh) or not np.array_equal(xyz[:, 2], ccx_mesh.coords_array[:, 2]):
            raise RuntimeError('mid/OML node-index correspondence differs')
    np.savez(out / 'mesh.npz', xyz=xyz, ccx_xyz=ccx_mesh.coords_array)
    result = {'reference': ref, 'placement': 'centred OML' if ref == 'oml' else
              'b: CCX on OML nodes with inward skin OFFSET; centred webs/TE_FLAT',
              'nodes': n, 'job': os.environ['SLURM_JOB_ID'], 'cases': {}}
    original = StaticLinearSolver._build_f_ext
    captured = {}

    def capture(self):
        f = original(self)
        captured['f'] = np.asarray(f).copy()
        return f

    cases = {'LC1': {'type': 'gravity', 'value': (0., -9.81, 0.)},
             'LC2': {'type': 'tip', 'value': (0., 1e6, 0.)},
             'LC4': {'type': 'distributed', 'value': (0., 5e3, 0.), 'n_slices': 60}}
    for name, load in cases.items():
        StaticLinearSolver._build_f_ext = capture
        try:
            u = run_shell_case(mesh, props, load).reshape(n, 6)
        finally:
            StaticLinearSolver._build_f_ext = original
        f = captured['f'].reshape(n, 6)
        if np.any(f[:, 3:]):
            raise RuntimeError('writer does not support nodal moments')
        np.savez(out / f'{name}-ae.npz', f=f, u=u)
        field = {node.id: f[i, :3].tolist() for i, node in enumerate(ccx_mesh.nodes)
                 if np.any(f[i, :3])}
        deck = out / f'{name}.inp'
        write_ccx_mesh(ccx_mesh, str(deck), properties=ccx_props, boundary_nodeset='RootNodes',
                       solver_type='LinearStatic', quadratic=True,
                       span_direction=(0., 0., 1.), load_field=field)
        # Suppress expensive unrequested stress output, preserve 2D displacement output.
        text = deck.read_text().replace('*EL FILE, OUTPUT=2D\nS, E\n', '')
        deck.write_text(text)
        if ref == 'mid':
            result['offset_sections'] = offset_sections(deck, ccx_mesh, ccx_props)
        fd = read_deck_force(deck, n)
        F = f[:, :3].sum(axis=0)
        Fd = fd[:, :3].sum(axis=0)
        err = float(np.linalg.norm(Fd-F)/np.linalg.norm(F))
        if err > 1e-6:
            raise RuntimeError(f'force transfer error {err}')
        row = {'ae': metric(xyz, f, u), 'Fae_N': F.tolist(), 'Fdeck_N': Fd.tolist(),
               'Mroot_Nm': np.cross(xyz, f[:, :3]).sum(axis=0).tolist(),
               'load_transfer_relative': err}
        result['cases'][name] = row
        save(out / 'result.json', result)
        p = run_ccx(deck, CCX)
        (out / f'{name}.stdout').write_text(p.stdout + p.stderr)
        if p.returncode:
            raise RuntimeError(f'CCX {name} rc={p.returncode}: {p.stdout[-2000:]}')
        disp = parse_frd_disp(deck.with_suffix('.frd'), list(range(1, n+1)))
        if len(disp) != n:
            raise RuntimeError(f'FRD matched {len(disp)}/{n}')
        uc = np.zeros((n, 6))
        for k, v in disp.items():
            uc[k-1, :3] = v
        row['ccx'] = metric(xyz, fd, uc)
        row['gap_pct_ae_over_ccx'] = {k: 100*(row['ae'][k]/row['ccx'][k]-1)
                                            for k in row['ae']}
        save(out / 'result.json', result)
        print(name, row, flush=True)

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily
    from run_s4_rotating_modal import assemble_rotating
    asm = MeshAssembler(mesh, {'elements': {'element_family': ElementFamily.SHELL,
                           'span_direction': (0., 0., 1.), 'properties': props}})
    freqs, modes, free = assemble_rotating(asm, mesh, 0., 8)
    full = np.zeros((8, n*6))
    full[:, free] = modes
    result['ae_modes'] = mode_rows(mesh, freqs, full.reshape(8, n, 6))
    save(out / 'result.json', result)
    deck = out / 'modal.inp'
    write_ccx_mesh(ccx_mesh, str(deck), properties=ccx_props, boundary_nodeset='RootNodes',
                   solver_type='Modal', num_modes=8, quadratic=True,
                   span_direction=(0., 0., 1.))
    deck.write_text(deck.read_text().replace('*EL FILE, OUTPUT=2D\nS, E\n', ''))
    if ref == 'mid':
        offset_sections(deck, ccx_mesh, ccx_props)
    p = run_ccx(deck, CCX)
    (out / 'modal.stdout').write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(f'CCX modal rc={p.returncode}: {p.stdout[-2000:]}')
    freqs = parse_ccx_frequencies(deck.with_suffix('.dat'), 8)
    modes = frd_modes(deck.with_suffix('.frd'), n)
    if len(modes) != len(freqs):
        raise RuntimeError('frequency/shape count mismatch')
    result['ccx_modes'] = mode_rows(mesh, freqs, modes)
    save(out / 'result.json', result)
    print(json.dumps(result, indent=2), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('reference', choices=['oml', 'mid', 'probe'])
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(REPO / 'logs/check/t10'):
        ap.error('--out must be under logs/check/t10')
    out.mkdir(parents=True, exist_ok=True)
    if args.reference == 'probe':
        probe(out)
    else:
        blade(args.reference, out)


if __name__ == '__main__':
    main()
