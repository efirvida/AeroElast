#!/usr/bin/env python3
"""T6/S6 independent S8R centrifugal shell arbiter; NOT an accuracy gate.

scripts/aeroenv.sh python -B tools/t6_ccx_centrifugal_buckling.py {oml,offset} --out logs/check/tipmodes/t6-ccx-REF
Blade runs require SLURM. Always es=0.25, official YAML, OML mesh; `offset`
reuses the T10 inward skin OFFSET (centred webs/TE_FLAT), not inset nodes.

The static step is linear (no NLGEOM). CENTRIF takes omega**2, not rpm or
omega: centre (0,0,0), axis +y, exactly as diagnose_rotating_prestress.py.
There is deliberately NO added hub radius: the diagnostic loads the generated
blade coordinates directly. Static global integration-point Szz is averaged
within each ply and thickness-integrated, then area-integrated over windows.
The reference is the diagnostic's discrete outboard nodal centrifugal load.
L8(b): any window outside 2% invalidates this arbiter BEFORE mode reporting.

BUCKLE without PERTURBATION uses the applied CENTRIF as the reference load,
not an incremental load atop the static preload (tiny smoke checks this).
FREQUENCY uses *STEP,PERTURBATION to retain the preceding static prestress.
CCX also linearises on the statically displaced geometry in that step; this
is an independent shell reference, not identical to our undeformed K+KG.
L8(c): physical only if lambda<=1.5, tip TE peak >=0.9 AND local Hz<3.8;
artifact if lambda>=3 and no such local mode; otherwise unresolved. No
frequency error tolerance or automatic physical/artifact verdict is imposed.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'tools'))
from run_s4_rotating_modal import AXIS, BEAM_FRACTION_MIN, beam_projection, omega_from_rpm  # noqa: E402  (after the sys.path insert)
from t10_ccx_arbiter import CCX, frd_modes, offset_sections, save  # noqa: E402

WINDOWS = ((0., .15), (.15, .3), (.3, .5), (.5, .7), (.7, .8), (.8, .9),
           (.9, .94), (.94, .97), (.97, 1.000001))


def centrifugal_card(rpm=7.56):
    """CCX CENTRIF magnitude is squared angular speed in rad^2/s^2."""
    return f'*DLOAD\nT6ALL,CENTRIF,{omega_from_rpm(rpm)**2:.16g},0,0,0,0,1,0\n'


def static_step():
    return ('*STEP\n*STATIC\n' + centrifugal_card() +
            '*EL PRINT,ELSET=T6ALL,GLOBAL=YES\nS\n' +
            '*NODE FILE,OUTPUT=2D\nU\n*END STEP\n')


def buckling_factors(path):
    values, active = [], False
    with path.open() as handle:
        for line in handle:
            if 'B U C K L I N G' in line:
                active = True
            elif active and 'E I G E N V A L U E' in line:
                break
            elif active and re.match(r'^\s*\d+\s+[+-]?[\d.]+[EDed][+-]?\d+\s*$', line):
                values.append(float(line.split()[1].replace('D', 'E')))
    if not values:
        raise RuntimeError(f'no buckling factors in {path}')
    return np.asarray(values)


def section_check(path, layers, area, cen, z, fz):
    """First static stress block only; no surface-stress/nodal extrapolation."""
    points = {}
    active = False
    with path.open() as handle:
        for line in handle:
            if 'stresses (elem, integ.pnt.' in line:
                if active:
                    break
                active = True
            elif active:
                fields = line.split()
                if len(fields) == 8 and fields[0].isdigit() and fields[1].isdigit():
                    eid, ip = int(fields[0]), int(fields[1])
                    points.setdefault(eid, {})[ip] = float(fields[4].replace('D', 'E'))
                elif points and line.strip():
                    break
    nzz = np.empty(len(area))
    for i, ts in enumerate(layers):
        p = points.get(i+1, {})
        if not p or len(p) % len(ts) or sorted(p) != list(range(1, len(p)+1)):
            raise RuntimeError(f'element {i+1}: invalid/missing integration points')
        per_ply = len(p) // len(ts)
        if per_ply not in (6, 8, 9):  # 6: triangles (S6 expanded to wedges)
            raise RuntimeError(f'element {i+1}: unexpected {per_ply} points/ply')
        stresses = np.asarray([p[j] for j in sorted(p)]).reshape(len(ts), per_ply)
        nzz[i] = np.dot(np.asarray(ts), stresses.mean(axis=1))
    z0, span = z.min(), np.ptp(z)
    zl = (cen[:, 2] - z0) / span
    rows = []
    for lo, hi in WINDOWS:
        a, b = z0 + lo*span, min(z0 + hi*span, z.max())
        # Integral_a^b sum_{z_node>zeta} fz_node dzeta, exactly (no centroid-bin aliasing).
        reference = float(np.dot(fz, np.clip(z-a, 0, b-a)))
        selected = (zl >= lo) & (zl < hi)
        fe = float(np.sum(nzz[selected] * area[selected]))
        ratio = fe/reference
        rows.append({'zL': [lo, hi], 'integral_ccx_Nm': fe,
                     'integral_ref_Nm': reference, 'ratio': ratio,
                     'within_2pct': bool(abs(ratio-1) <= .02)})
    return rows, nzz


def shape_rows(mesh, values, shapes, kind):
    xyz = mesh.coords_array
    names = {}
    for name in sorted(mesh.element_sets):
        for el in mesh.element_sets[name].elements:
            for nid in el.node_ids:
                names.setdefault(nid, set()).add(name)
    rows = []
    for value, shape in zip(values, shapes, strict=True):
        u = shape[:, :3]
        frac, _, residual = beam_projection(xyz, u)
        peak = int(np.argmax(np.linalg.norm(u, axis=1)))
        rpeak = int(np.argmax(np.linalg.norm(residual, axis=1)))
        zl = (xyz[peak, 2]-xyz[:, 2].min())/np.ptp(xyz[:, 2])
        rzl = (xyz[rpeak, 2]-xyz[:, 2].min())/np.ptp(xyz[:, 2])
        rows.append({kind: float(value), 'beam_fraction': frac,
                     'classification': 'local' if frac < BEAM_FRACTION_MIN else 'beam',
                     'peak_zL': float(zl), 'residual_peak_zL': float(rzl),
                     'peak_sets': sorted(names.get(mesh.nodes[peak].id, set()))})
    return rows


def run(deck):
    from tests.support.ccx_io import run_ccx
    p = run_ccx(deck, CCX)
    deck.with_suffix('.stdout').write_text(p.stdout + p.stderr)
    if p.returncode or '*ERROR' in p.stdout:
        raise RuntimeError(f'CCX {deck} rc={p.returncode}; see {deck.with_suffix(".stdout")}')


def blade(realisation, out, stage='all'):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('es=0.25 blade runs require sbatch')
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.core.mesh.io.writers import write_ccx_mesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=.25, skin_reference='oml',
                    yaml_file=str(REPO / 'tests/IEA-15-240-RWT.yaml'))
    mesh = gen.generate(renumber='rcm')
    props = build_rust_properties(gen.numad_mesh_data)
    deck = out / 'template.inp'
    write_ccx_mesh(mesh, str(deck), properties=props, boundary_nodeset='RootNodes',
                   solver_type='LinearStatic', quadratic=True, span_direction=(0., 0., 1.))
    offsets = offset_sections(deck, mesh, props) if realisation == 'offset' else {}
    # Keep the shared writer's material/orientation/sets/expanded-root invariants.
    model = deck.read_text().split('*STEP\n', 1)[0]
    model += f'*ELSET,ELSET=T6ALL,GENERATE\n1,{len(mesh.elements)},1\n'
    static = static_step()
    (out / 'static.inp').write_text(model + static)
    (out / 'buckle.inp').write_text(model + static + '*STEP\n*BUCKLE\n20\n' +
                                      centrifugal_card() + '*EL PRINT,ELSET=T6ALL,FREQUENCY=0\nS\n' +
                                      '*NODE FILE,OUTPUT=2D\nU\n*END STEP\n')
    (out / 'frequency.inp').write_text(model + static + '*STEP,PERTURBATION\n*FREQUENCY\n15\n' +
                                         '*EL PRINT,ELSET=T6ALL,FREQUENCY=0\nS\n' +
                                         '*NODE FILE,OUTPUT=2D\nU\n*END STEP\n')
    asm = MeshAssembler(mesh, {'elements': {'element_family': ElementFamily.SHELL,
                        'span_direction': (0., 0., 1.), 'properties': props}})
    xyz = mesh.coords_array
    fcf = asm.assemble_centrifugal_load(omega_from_rpm(7.56), AXIS, np.zeros(3))
    layers_by_id = {}
    for sec in gen.numad_mesh_data['sections']:
        for el in mesh.element_sets[sec['elementSet']].elements:
            layers_by_id.setdefault(el.id, [float(t) for _, t, _ in sec['layup']])
    cen, area, layers = [], [], []
    for el in mesh.elements:
        p = np.asarray([n.coords for n in el.nodes])
        cen.append(p.mean(axis=0))
        area.append(.5*np.linalg.norm(np.cross(p[2]-p[0], p[3]-p[1])) if len(p) == 4 else
                    .5*np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0])))
        layers.append(layers_by_id[el.id])
    area, cen = np.asarray(area), np.asarray(cen)
    result = {'realisation': realisation, 'es': .25, 'rpm': 7.56, 'axis': AXIS.tolist(),
              'center_m': [0., 0., 0.], 'added_hub_offset_m': 0.,
              'job_id': os.environ['SLURM_JOB_ID'], 'offset_sections': offsets}
    np.savez_compressed(out / 'mesh-load.npz', xyz=xyz, cen=cen, area=area, f_cf=fcf)
    if stage in ('buckle', 'frequency'):
        # Split stages (sequana_cpu_dev, 20 min each): an eigen stage needs a prior valid static stage.
        result = json.loads((out / 'result.json').read_text())
        if not result.get('prestress_valid'):
            raise RuntimeError('run --stage static first; L8(b) prestress check missing or failed')
        result['job_id_' + stage] = os.environ['SLURM_JOB_ID']
        eigen_stage(out, stage, mesh, xyz, result)
        return
    run(out / 'static.inp')
    rows, nzz = section_check(out / 'static.dat', layers, area, cen, xyz[:, 2], fcf[2::6])
    result['prestress_windows'] = rows
    result['prestress_valid'] = all(r['within_2pct'] for r in rows)
    save(out / 'result.json', result)
    np.savez_compressed(out / 'prestress.npz', nzz=nzz, cen=cen, area=area)
    # Report this BEFORE either eigen analysis; never classify an invalid CENTRIF setup.
    print('PRESTRESS CHECK (L8b):', json.dumps(rows), flush=True)
    if not result['prestress_valid']:
        raise RuntimeError('L8(b) section tension outside 2%; arbiter invalid, eigen runs withheld')
    if stage == 'static':
        return
    for kind in ('buckle', 'frequency'):
        eigen_stage(out, kind, mesh, xyz, result)


def eigen_stage(out, kind, mesh, xyz, result):
    from tests.support.ccx_io import parse_ccx_frequencies
    deck = out / f'{kind}.inp'
    run(deck)
    values = (buckling_factors(deck.with_suffix('.dat')) if kind == 'buckle' else
              parse_ccx_frequencies(deck.with_suffix('.dat'), 15))
    # OUTPUT=2D must be set already in the static step (CCX 2.20).
    # The leading DISP blocks are the preload (BUCKLE also writes a time=0 reference block);
    # the eigenmodes are always the last len(values) blocks.
    blocks = frd_modes(deck.with_suffix('.frd'), len(xyz))
    if len(blocks) < len(values) + 1:
        raise RuntimeError(f'{kind}: expected preload block(s) + {len(values)} modes, got {len(blocks)}')
    shapes = blocks[-len(values):]
    if len(shapes) != len(values):
        raise RuntimeError(f'{kind}: {len(shapes)} shapes for {len(values)} eigenvalues')
    key = 'lambda' if kind == 'buckle' else 'Hz'
    result[kind] = shape_rows(mesh, values, shapes, key)
    np.savez_compressed(out / f'{kind}.npz', values=values, modes=np.asarray(shapes), xyz=xyz)
    save(out / 'result.json', result)
    print(kind, json.dumps(result[kind]), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('realisation', choices=['oml', 'offset'])
    ap.add_argument('--stage', choices=['all', 'static', 'buckle', 'frequency'], default='all',
                    help='run one stage per job (each fits the 20 min dev limit); eigen stages need static')
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    out = args.out.resolve()
    if not out.is_relative_to(REPO / 'logs/check/tipmodes'):
        ap.error('--out must be under logs/check/tipmodes')
    out.mkdir(parents=True, exist_ok=True)
    blade(args.realisation, out, args.stage)


if __name__ == '__main__':
    main()
