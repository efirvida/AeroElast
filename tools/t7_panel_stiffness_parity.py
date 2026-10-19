#!/usr/bin/env python3
"""T7(iv)/S7: parked MITC4 versus S8R panel stiffness, NOT an accuracy gate.

scripts/aeroenv.sh python -B tools/t7_panel_stiffness_parity.py --reference-only
scripts/aeroenv.sh python -B tools/t7_panel_stiffness_parity.py

Panel runs require SLURM; official IEA-15, es=0.25, 09_46_LP_TE_PANEL.
No K_G, rotation, prestress, pressure quadrature, or full-blade solve.
OML coordinates are used as the common laminate CENTRE, OFFSET=0, not as
an inward-stacked skin. Both codes project span_direction=(0,0,1) into each
element plane; CCX uses the shared writer's 0.1-degree angle buckets.
Clamp the complete inboard cut, all six DOFs, including CCX midsides.
Apply a 1 N resultant in the area-averaged normal direction, equally over
FREE CORNER nodes only. This is a discrete load, NOT physical pressure.
S8R midsides carry zero applied force. Work is 0.5*f.u over these same nodes.
Use the writer's serialized CLOAD values in BOTH solves. Re-emit only the
writer-generated mesh coordinate fields at full precision, retaining its
labels/connectivity and evaluating midsides on the original straight edges.

Before any panel solve, validate the CCX path against a homogeneous axial
cantilever (nu=0, L=1, b=.1, h=.001, E=210 GPa, F=1 N):
(A) work versus F^2 L/(2 E b h), using the existing 2% membrane analytical
band in tests/validation/parity/test_composite_beam_parity.py;
(B) repeat at doubled mesh density, report work change and check the same
analytical band. Also check load doubling gives four times the work (1e-5
numerical/FRD precision check, NOT a physical parity band).
This validates axial reference/load/BC plumbing, not curved-panel bending.
Panel ratios have NO accepted physical tolerance and never print parity PASS.
Artifacts include decks, stdout, geometry, connectivity, discrete f/u and JSON.
Default output must be absent or empty; use --out for a fresh managed directory.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
MANAGED = REPO / 'logs/check/tipmodes/t7-iv'
SPAN = (0., 0., 1.)
PANEL = '09_46_LP_TE_PANEL'


def boundary_edges(mesh):
    """Topological boundary; reject non-quads and nonmanifold edges."""
    counts = Counter()
    for el in mesh.elements:
        if el.element_type.name != 'quad':
            raise ValueError('patch must contain only linear quads')
        ids = el.node_ids
        counts.update(tuple(sorted((ids[i], ids[(i+1) % 4]))) for i in range(4))
    if any(n > 2 for n in counts.values()):
        raise ValueError('nonmanifold patch')
    return [edge for edge, n in counts.items() if n == 1]


def prepare_patch(source, name, span=SPAN):
    """Use public extraction; preserve connectivity, property sets and winding."""
    from aeroelast.core.mesh.entities import NodeSet

    mesh = source.extract_submesh(name)
    mesh.renumber_mesh('simple')
    boundary = boundary_edges(mesh)
    # Edge-connectedness is required, not merely touching at a corner.
    pending = set(range(len(mesh.elements)))
    reached = {pending.pop()}
    edges = [{tuple(sorted((e.node_ids[i], e.node_ids[(i+1) % 4])))
              for i in range(4)} for e in mesh.elements]
    while pending:
        new = {i for i in pending if any(edges[i] & edges[j] for j in reached)}
        if not new:
            raise ValueError('patch must be edge-connected')
        pending -= new
        reached |= new
    s = mesh.coords_array @ np.asarray(span)
    tol = max(float(np.ptp(s)), 1.) * 1e-8
    fixed_ids = {n for edge in boundary for n in edge
                 if all(abs(s[mesh.node_id_to_index[k]] - s.min()) <= tol for k in edge)}
    if len(fixed_ids) < 2:
        raise ValueError('no complete inboard boundary edge at minimum span coordinate')
    mesh.add_node_set(NodeSet('T7FIXED', [mesh.node_map[n] for n in fixed_ids]))
    fixed = sorted(mesh.node_id_to_index[n] for n in fixed_ids)
    free = np.setdiff1d(np.arange(len(mesh.nodes)), fixed)
    normals = [np.cross(e.node_coords[2]-e.node_coords[0],
                        e.node_coords[3]-e.node_coords[1]) for e in mesh.elements]
    normal = np.sum(normals, axis=0)
    if not len(free) or np.linalg.norm(normal) < 1e-14:
        raise ValueError('patch has no free nodes or usable mean normal')
    if any(np.dot(n, normal) <= 0 for n in normals):
        raise ValueError('patch winding is inconsistent')
    force = np.zeros((len(mesh.nodes), 6))
    force[free, :3] = normal / np.linalg.norm(normal) / len(free)
    return mesh, fixed, force.ravel()


def write_deck(mesh, props, force, deck, span=SPAN):
    """Reuse writer's layup, orientation buckets, nodal labels and expanded BC."""
    from aeroelast.core.mesh.io.writers import write_ccx_mesh

    field = {n.id: force.reshape(-1, 6)[i, :3].tolist() for i, n in enumerate(mesh.nodes)}
    write_ccx_mesh(mesh, str(deck), properties=props, boundary_nodeset='T7FIXED',
                   solver_type='LinearStatic', quadratic=True, span_direction=span,
                   load_field=field)
    # The shared writer prints 7 significant digits. Keep its topology and angle
    # buckets, but avoid rounding a 100 m span coordinate before a .1 m patch solve.
    path = deck.with_suffix('.msh')
    lines = path.read_text().splitlines()
    xyz = {i+1: n.coords for i, n in enumerate(mesh.nodes)}
    in_elements = False
    for line in lines:
        if line.startswith('*'):
            in_elements = line.upper().startswith('*ELEMENT')
        elif in_elements and line.strip():
            labels = [int(v) for v in line.split(',')]
            corners, midsides = labels[1:5], labels[5:9]
            for i, mid in enumerate(midsides):
                xyz[mid] = (xyz[corners[i]] + xyz[corners[(i+1) % 4]]) / 2
    in_nodes = False
    for i, line in enumerate(lines):
        if line.startswith('*'):
            in_nodes = line.upper().startswith('*NODE')
        elif in_nodes and line.strip():
            label = int(line.split(',')[0])
            lines[i] = f'{label},' + ','.join(f'{v:.17g}' for v in xyz[label])
    path.write_text('\n'.join(lines)+'\n')
    # Read the public deck, never duplicate the writer's force-formatting rule.
    serialized = np.zeros_like(force).reshape(-1, 6)
    active = False
    for line in deck.read_text().splitlines():
        if line.startswith('*'):
            active = line.upper() == '*CLOAD'
        elif active and line.strip():
            node, dof, value = line.split(',')
            serialized[int(node)-1, int(dof)-1] += float(value)
    return serialized.ravel()


def solve_mitc4(mesh, props, fixed, force, span=SPAN):
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import spsolve

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily

    asm = MeshAssembler(mesh, {'elements': {'element_family': ElementFamily.SHELL,
                        'span_direction': span, 'properties': props}})
    mat = asm.assemble_stiffness_matrix()
    if mat.comm.size != 1:
        raise ValueError('run serial, not with mpirun')
    indptr, indices, values = mat.getValuesCSR()
    k = csr_matrix((values.copy(), indices.copy(), indptr.copy()), shape=mat.getSize())
    mat.destroy()
    constrained = np.asarray([6*i+d for i in fixed for d in range(6)])
    free = np.setdiff1d(np.arange(len(force)), constrained)
    u = np.zeros_like(force)
    u[free] = spsolve(k[free][:, free], force[free])
    residual = k @ u - force
    work = float(.5 * force @ u)
    result = {'work_J': work,
              'free_residual': float(np.linalg.norm(residual[free]) / np.linalg.norm(force)),
              'energy_relative_error': float(abs(.5*u @ (k @ u)-work) / abs(work)),
              'force_balance_relative_error': float(np.linalg.norm(
                  residual.reshape(-1, 6)[fixed, :3].sum(axis=0) +
                  force.reshape(-1, 6)[:, :3].sum(axis=0)) / np.linalg.norm(force))}
    if not np.all(np.isfinite(u)) or work <= 0 or result['free_residual'] > 1e-7:
        raise RuntimeError(f'invalid MITC4 static solve: {result}')
    return result, u


def solve_ccx(mesh, props, force, deck, ccx, span=SPAN):
    from tests.support.ccx_io import parse_frd_disp, run_ccx

    serialized = write_deck(mesh, props, force, deck, span)
    if not np.array_equal(force, serialized):
        raise ValueError('force must first be canonicalized through write_deck')
    run = run_ccx(deck, ccx)
    deck.with_suffix('.stdout').write_text(run.stdout + run.stderr)
    if run.returncode or '*ERROR' in run.stdout.upper():
        raise RuntimeError(f'CCX failed rc={run.returncode}; see {deck.with_suffix(".stdout")}')
    labels = list(range(1, len(mesh.nodes)+1))
    disp = parse_frd_disp(deck.with_suffix('.frd'), labels)
    if set(disp) != set(labels):
        raise RuntimeError('CCX missing corner displacement(s)')
    u = np.asarray([disp[i] for i in labels])
    work = float(.5 * np.sum(force.reshape(-1, 6)[:, :3] * u))
    if not np.isfinite(work) or work <= 0:
        raise RuntimeError('CCX nonpositive or nonfinite load work')
    return work, u


def reference(out, ccx):
    """Two mesh densities and a load-scaling check BEFORE panel comparison."""
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
    from aeroelast.core.mesh.model import MeshModel
    from aeroelast.models.blade.model import build_rust_properties

    # Two identical plies exercise the COMPOSITE section path used on the blade.
    props = build_rust_properties({'materials': [{'name': 'reference', 'density': 7800.,
        'elastic': {'E': 210e9, 'nu': 0.}}], 'sections': [{'elementSet': 'ref',
        'layup': [['reference', .0005, 0.], ['reference', .0005, 0.]]}]})
    expected = 1. / (2 * 210e9 * .1 * .001)
    rows = []
    for nx, ny in ((4, 16), (8, 32)):
        mesh = MeshModel()
        grid = {}
        for j in range(ny+1):
            for i in range(nx+1):
                node = Node([.1*i/nx, j/ny, 0.])
                mesh.add_node(node)
                grid[i, j] = node
        for j in range(ny):
            for i in range(nx):
                mesh.add_element(MeshElement([grid[i, j], grid[i+1, j],
                    grid[i+1, j+1], grid[i, j+1]], ElementType.quad))
        mesh.add_element_set(ElementSet('ref', mesh.elements))
        mesh.add_node_set(NodeSet('T7FIXED', [grid[i, 0] for i in range(nx+1)]))
        mesh.renumber_mesh('simple')
        force = np.zeros((len(mesh.nodes), 6))
        for i in range(nx+1):
            force[mesh.node_id_to_index[grid[i, ny].id], 1] = (.5 if i in (0, nx) else 1.)/nx
        force = force.ravel()
        work, _ = solve_ccx(mesh, props, force, out / f'reference-{nx}.inp', ccx, (0, 1, 0))
        scaled, _ = solve_ccx(mesh, props, 2*force, out / f'doubled-{nx}.inp', ccx, (0, 1, 0))
        rows.append({'nx': nx, 'ny': ny, 'work_J': work, 'analytic_work_J': expected,
                     'relative_error': work/expected-1, 'doubling_work_ratio': scaled/work})
    result = {'checks': rows, 'mesh_work_change': rows[1]['work_J']/rows[0]['work_J']-1,
              'valid': all(abs(r['relative_error']) <= .02 and
                           abs(r['doubling_work_ratio']/4-1) <= 1e-5 for r in rows)}
    (out / 'reference.json').write_text(json.dumps(result, indent=2)+'\n')
    print('REFERENCE', json.dumps(result), flush=True)
    if not result['valid']:
        raise RuntimeError('reference invalid; panel comparison withheld')
    return result


def section_json_value(value):
    """Encode NumPy section metadata without masking unsupported report types."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f'Object of type {type(value).__name__} is not JSON serializable')


def panel(out, ccx, checked):
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=.25, skin_reference='oml',
                    yaml_file=str(REPO / 'tests/IEA-15-240-RWT.yaml'))
    source = gen.generate(renumber='rcm')
    mesh, fixed, force = prepare_patch(source, PANEL)
    props = {PANEL: build_rust_properties(gen.numad_mesh_data)[PANEL]}
    # Keep the exact section source alongside the exported material/orientation cards.
    sections = [s for s in gen.numad_mesh_data['sections'] if s['elementSet'] == PANEL]
    requested_force = force.copy()
    force = write_deck(mesh, props, force, out / 'panel.inp')
    ae, u = solve_mitc4(mesh, props, fixed, force)
    work, ccx_u = solve_ccx(mesh, props, force, out / 'panel.inp', ccx)
    np.savez_compressed(out / 'panel.npz', xyz=mesh.coords_array,
                        connectivity=np.array([e.node_ids for e in mesh.elements]),
                        fixed_indices=fixed, force=force, mitc4_u=u, ccx_u=ccx_u)
    result = {'reference': checked, 'panel': PANEL, 'es_m': .25, 'skin_reference': 'oml',
              'laminate_reference': 'centered on OML coordinates; OFFSET=0 in both codes',
              'span_direction': SPAN, 'sections': sections, 'nodes': len(mesh.nodes),
              'elements': len(mesh.elements), 'boundary': 'minimum-span complete cut, all six DOFs',
              'load': '1 N mean-normal resultant, equal on free corners, zero on CCX midsides',
              'load_serialization_relative_change': float(np.linalg.norm(force-requested_force) /
                                                          np.linalg.norm(requested_force)),
              'resultant_N': force.reshape(-1, 6)[:, :3].sum(axis=0).tolist(),
              'mitc4': ae, 'ccx_work_J': work, 'mitc4_over_ccx_work': ae['work_J']/work,
              'physical_parity_band': None, 'verdict': 'diagnostic only'}
    (out / 'result.json').write_text(json.dumps(result, indent=2, default=section_json_value)+'\n')
    print('PANEL', json.dumps(result, default=section_json_value), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference-only', action='store_true', help='validate CCX plumbing without blade meshing')
    ap.add_argument('--out', type=Path, default=MANAGED / 'run', help='fresh managed output directory')
    ap.add_argument('--ccx', default=os.environ.get('CCX_BIN', shutil.which('ccx') or
                    '/scratch/leahk/eduardo.donestevez/venv/bin/ccx'), help='CalculiX executable')
    args = ap.parse_args(argv)
    out = args.out.resolve()
    if not out.is_relative_to(MANAGED.resolve()):
        ap.error('--out must be under logs/check/tipmodes/t7-iv')
    if not args.reference_only and not os.environ.get('SLURM_JOB_ID'):
        ap.error('panel comparison requires SLURM; use --reference-only locally')
    if out.exists() and any(out.iterdir()):
        ap.error('--out must be absent or empty')
    if not shutil.which(args.ccx):
        ap.error('CalculiX executable not found')
    out.mkdir(parents=True, exist_ok=True)
    checked = reference(out, args.ccx)
    if not args.reference_only:
        panel(out, args.ccx, checked)


if __name__ == '__main__':
    main()
