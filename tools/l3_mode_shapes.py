#!/usr/bin/env python
"""Mode shapes of the IEA-15 blade: shell (es 0.25, mid skin) vs BeamDyn p20, MAC per mode (l3-closure T8, S6).

Three steps, each writes under ``logs/check/l3/t8``:

  bd      BeamDyn p20 eigenvectors from the existing linearizations (no OpenFAST run), parked and rated; gate: the
          frequencies reproduce ``logs/check/of_bd/<case>/modes.json`` to 1e-6.
  shell   the S-4 modal problem (``tests/test_iea15mw_s4_rotating_modal.py::_build_assembly`` and
          ``tools/run_s4_rotating_modal.assemble_rotating``) with the eigenvectors saved (``shell-{parked,rated}.npz``).
  compare section-mean translations and twist at the BeamDyn node heights, MAC, Hungarian pairing, ``RESULT.md``.

BeamDyn state -> node kinematics (``bd_modes_from_lin``). The continuous states of ``case.*.lin`` are, per blade, the
finite-element node displacements (derivative order 2: ``BD_1 finite element node N ... translational displacement in
X, m`` and ``rotational displacement in X, rad``, nodes 2..21; node 1 is the clamped root) followed by their time
derivatives (``BD_1 First time derivative of ...``). The MBC3 collective block keeps the first state of each blade
triplet (``RotTripletIndicesStates2[:, 0]``), so ``vectors[:len(q), j]`` of the eigenvector of ``[q, qdot]`` is the
displacement part: translations of the node reference point and the rotation of the node relative to its reference
orientation, both in the blade-root (BeamDyn global) axes. The rotation states are Wiener-Milenkovic parameters
``c = 4 tan(phi/4) n``, equal to the rotation vector to first order, which is what a linear mode contains. Rotation
about BeamDyn z is the twist about the span axis (z is common to both frames).

Node heights: the p20 member has GLL nodes in the normalised coordinate of the member; ``z = L (1 + xi) / 2`` with
``L = 117 m`` reproduces the nine tip nodes of the p60 BeamDyn run (initial VTK, ``bd_static_reference.py``) to 1e-4 m
(``test_gll_node_heights_match_the_p60_beamdyn_nodes``).

Frame: BeamDyn (x, y, z) -> shell (-y, x, z) (``tools/bd_static_reference.py``): shell x is the chord (edge), shell y the
flap direction.

Shell station fields (``station_fields``): per node ring (``tools/s2_ring_metric_blade.ring_profile``: nodes at one
z level, arc-length weights) the arc-weighted mean of (u_x, u_y) and the twist theta_z from ``rigid_rotation_fit`` on
the ring nodes (``fit='all'``) or on its box nodes (``fit='box'``, ``box_node_mask``; NaN below 6 nodes), linear in z
between the two bracketing rings.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from tests.support.modal_pairing import mac_matrix, pair_modes_by_mac  # noqa: E402

OUT = ROOT / "logs/check/l3/t8"
BD_LOGS = ROOT / "logs/check/of_bd"
BLADE_LENGTH = 117.0
RPM_RATED = 7.56
N_MODES = 8
BD_CASES = {"parked": "bp20_0", "rated": "bp20_7.56"}


def gll_nodes(order: int) -> np.ndarray:
    """The ``order + 1`` Gauss-Lobatto-Legendre nodes on [-1, 1]."""
    from numpy.polynomial import legendre as leg

    c = np.zeros(order + 1)
    c[-1] = 1.0
    return np.concatenate([[-1.0], np.sort(leg.legroots(leg.legder(c))), [1.0]])


def bd_node_z(order: int, length: float = BLADE_LENGTH) -> np.ndarray:
    """Initial z of BeamDyn nodes 1..order+1 of a single member of ``order``: GLL points over the blade length."""
    return length * (1.0 + gll_nodes(order)) / 2.0


def real_shape(v: np.ndarray) -> np.ndarray:
    """Rotate a complex mode to its dominant real phase, take the real part, normalise, fix the sign.

    ``phi = angle(sum v^2) / 2`` maximises ``sum Re(v e^{-i phi})^2``; the sign is fixed so the largest component is
    positive. The result is the same for ``v``, ``-v`` and ``v e^{i alpha}`` (up to rounding).
    """
    v = np.asarray(v, dtype=complex)
    phi = 0.5 * np.angle(np.sum(v * v))
    r = (v * np.exp(-1j * phi)).real
    r = r / max(np.linalg.norm(r), 1e-300)
    return r * (1.0 if r[np.argmax(np.abs(r))] >= 0 else -1.0)


def pair_shell_to_bd(shell: np.ndarray, bd: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hungarian MAC pairing; ``pair[j]`` is the shell row of BeamDyn mode ``j``; ``mac`` is (n_shell, n_bd)."""
    return pair_modes_by_mac(np.asarray(shell), np.asarray(bd))


_STATE = re.compile(r"finite element node (\d+) .*element order = (\d+)\) (translational|rotational) displacement in ([XYZ])")


def bd_modes_from_lin(lin_files, n_modes: int = N_MODES) -> dict:
    """Frequencies and node translations / rotations of the BeamDyn collective modes (as ``post()`` in
    ``tools/of_bd_blade_modal.py``: MBC3, blade block of ``AvgA``, ``eig``, ``imag > 0`` sorted by ``|ev|``).

    Returns ``freqs`` (Hz, ``|ev| / 2 pi``), ``trans`` and ``rot`` (n_modes, n_nodes, 3) complex in BeamDyn axes,
    ``z`` (n_nodes,), ``nodes`` (BeamDyn node numbers).
    """
    from openfast_toolbox.linearization.mbc import fx_mbc3

    mbc, data = fx_mbc3([str(p) for p in lin_files], verbose=False)
    if not mbc["performedTransformation"]:
        raise RuntimeError("MBC3 transformation not performed")
    q = np.asarray(data["RotTripletIndicesStates2"], dtype=int)[:, 0]
    if data["ndof1"] or 3 * len(q) != data["ndof2"]:
        raise RuntimeError("Expected blade-only second-order rotating states")
    idx = np.concatenate([q, q + data["ndof2"]])
    ev, vec = np.linalg.eig(mbc["AvgA"][np.ix_(idx, idx)])
    desc = [data["DescStates"][k] for k in q]
    parsed = []
    for d in desc:
        m = _STATE.search(d)
        if m is None or not d.lstrip().startswith("BD_1 finite"):
            raise RuntimeError(f"state is not a BD_1 node displacement: {d!r}")
        parsed.append((int(m[1]), int(m[2]), m[3], "XYZ".index(m[4])))
    nodes = sorted({p[0] for p in parsed})
    order = parsed[0][1]
    col = {n: i for i, n in enumerate(nodes)}
    sel = [j for j in np.argsort(np.abs(ev)) if ev[j].imag > 0][:n_modes]
    trans = np.zeros((len(sel), len(nodes), 3), dtype=complex)
    rot = np.zeros_like(trans)
    for m, j in enumerate(sel):
        for k, (node, _o, kind, comp) in enumerate(parsed):
            (trans if kind == "translational" else rot)[m, col[node], comp] = vec[k, j]
    return {
        "freqs": np.array([abs(ev[j]) / (2 * np.pi) for j in sel]),
        "trans": trans,
        "rot": rot,
        "z": bd_node_z(order)[np.array(nodes) - 1],
        "nodes": np.array(nodes),
    }


def bd_to_shell_axes(trans: np.ndarray, rot: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """BeamDyn (x, y, z) -> shell (-y, x, z) for translations and rotation vectors."""
    swap = lambda a: np.stack([-a[..., 1], a[..., 0], a[..., 2]], axis=-1)  # noqa: E731
    return swap(trans), swap(rot)


def station_fields(coords, u, ring_id, w, box, z_st, fit: str = "all"):
    """(u_x, u_y, theta_z) of the rings at the stations ``z_st``: arc-weighted mean translation and rigid-fit twist.

    ``ring_id`` (n,) is the ring index of each node (-1: none), ``w`` the arc weights, ``box`` the box-node mask,
    ``u`` (n, >=3) the nodal translations. Linear interpolation in z between the bracketing rings over the rings with
    a finite value; a station outside them is NaN.
    """
    from aeroelast.postprocess.sectional import station_rigid_rotation

    coords = np.asarray(coords)
    u = np.asarray(u)
    n_ring = int(ring_id.max()) + 1
    zr = np.full(n_ring, np.nan)
    mean = np.full((n_ring, 2), np.nan)
    theta = np.full(n_ring, np.nan)
    for k in range(n_ring):
        nodes = np.flatnonzero((ring_id == k) & (w > 0))
        if nodes.size == 0:
            continue
        zr[k] = coords[nodes, 2].mean()
        mean[k] = (u[nodes, :2] * w[nodes, None]).sum(axis=0) / w[nodes].sum()
        theta[k] = station_rigid_rotation(coords, nodes, u[:, :3], fit_nodes=fit, box_mask=box)[2]
    order = np.argsort(zr)
    zr, mean, theta = zr[order], mean[order], theta[order]

    def interp(y):
        ok = np.isfinite(y) & np.isfinite(zr)
        if not ok.any():
            return np.full(len(z_st), np.nan)
        out = np.interp(z_st, zr[ok], y[ok])
        out[(z_st < zr[ok].min() - 1e-9) | (z_st > zr[ok].max() + 1e-9)] = np.nan
        return out

    return interp(mean[:, 0]), interp(mean[:, 1]), interp(theta)


def chord_at(z, yaml_file: Path = ROOT / "tests/IEA-15-240-RWT.yaml") -> np.ndarray:
    """Local chord (m) of the official yaml at ``z`` (the yaml chord grid is the fraction z / 117)."""
    import yaml

    c = yaml.safe_load(yaml_file.read_text())["components"]["blade"]["outer_shape_bem"]["chord"]
    return np.interp(np.asarray(z) / BLADE_LENGTH, c["grid"], c["values"])


def vectors(ux, uy, th, chord):
    """(bending, twist, combined) from station arrays, each (n_modes, n_st); rotation scaled by chord / 2."""
    bend = np.concatenate([ux, uy], axis=1)
    return bend, th, np.concatenate([ux, uy, th * chord / 2.0], axis=1)


def mac_report(shell_vec, bd_vec):
    """Pair on the combined vector; MAC of the three vectors for the pair."""
    pair, _ = pair_shell_to_bd(shell_vec[2], bd_vec[2])
    macs = [mac_matrix(s, b) for s, b in zip(shell_vec, bd_vec, strict=True)]
    return pair, macs


def family(ux, uy, th, chord):
    """Dominant kinematic family of one real mode and the energy fractions (edge x, flap y, twist theta c / 2)."""
    e = np.array([np.sum(ux**2), np.sum(uy**2), np.sum((th * chord / 2) ** 2)])
    e = e / e.sum()
    return ("E", "F", "T")[int(np.argmax(e))], e


def family_labels(ux, uy, th, chord) -> list[str]:
    """``1F``, ``1E``, ``2F``, ... of the modes (rows of ``ux``, ``uy``, ``th``), counted per family in row order."""
    counts = {"F": 0, "E": 0, "T": 0}
    labels = []
    for m in range(len(ux)):
        fm, _e = family(ux[m], uy[m], th[m], chord)
        counts[fm] += 1
        labels.append(f"{counts[fm]}{fm}")
    return labels


def bd_station_modes(bd, n_modes: int = N_MODES):
    """Real station shapes of the first ``n_modes`` BeamDyn modes in shell axes: ``(z, chord, (u_x, u_y, theta))``.

    ``bd`` is the ``bd_modes_from_lin`` dict or its ``.npz``. Each complex mode is turned to the dominant real phase of
    its combined vector ``[u_x; u_y; theta c / 2]`` (one common phase for the three parts), then its real part is kept.
    """
    z = np.asarray(bd["z"])
    chord = chord_at(z)
    t, r = bd_to_shell_axes(np.asarray(bd["trans"]), np.asarray(bd["rot"]))
    ux, uy, th = (np.zeros((n_modes, len(z))) for _ in range(3))
    for m in range(n_modes):
        comb = np.concatenate([t[m, :, 0], t[m, :, 1], r[m, :, 2] * chord / 2])
        ph = np.exp(-0.5j * np.angle(np.sum(comb * comb)))
        ux[m], uy[m], th[m] = (t[m, :, 0] * ph).real, (t[m, :, 1] * ph).real, (r[m, :, 2] * ph).real
    return z, chord, (ux, uy, th)


def shell_station_modes(coords, modes, ring_id, w, box, z_st, fit: str = "all"):
    """``(u_x, u_y, theta)`` of every shell mode (``modes`` (n_modes, n_nodes, >=3)) at the stations, by ``station_fields``.

    A NaN station (box fit with fewer than 6 nodes, or outside the rings) counts zero; the number of modes that had one is
    returned as the fourth item so the caller reports it.
    """
    n = len(modes)
    ux, uy, th = (np.zeros((n, len(z_st))) for _ in range(3))
    for m in range(n):
        ux[m], uy[m], th[m] = station_fields(coords, modes[m], ring_id, w, box, z_st, fit=fit)
    bad = int(np.isnan(th).any(axis=1).sum())
    return ux, uy, np.nan_to_num(th), bad


# ---------------------------------------------------------------- steps


def step_bd():
    OUT.mkdir(parents=True, exist_ok=True)
    for key, case in BD_CASES.items():
        files = sorted((BD_LOGS / case).glob("case.*.lin"), key=lambda p: int(p.stem.split(".")[-1]))
        t0 = time.time()
        res = bd_modes_from_lin(files, n_modes=12)
        rec = [m["hz"] for m in json.loads((BD_LOGS / case / "modes.json").read_text())["modes"]][:12]
        rel = np.abs(res["freqs"] / np.array(rec) - 1)
        print(f"bd {key}: {len(files)} lin, max rel diff vs modes.json {rel.max():.2e} ({time.time() - t0:.1f} s)")
        if rel.max() > 1e-6:
            raise SystemExit(f"GATE FAILED: {key} frequencies differ from modes.json by {rel.max():.2e}")
        np.savez(OUT / f"bd-{key}.npz", **res)


def step_shell(key):
    from run_s4_rotating_modal import assemble_rotating, omega_from_rpm
    from s2_ring_metric_blade import ring_profile
    from tests.test_iea15mw_s4_rotating_modal import _build_assembly
    from aeroelast.postprocess.sectional import box_node_mask

    rpm = 0.0 if key == "parked" else RPM_RATED
    t0 = time.time()
    asm, mesh = _build_assembly(0.25)
    print(f"mesh {mesh.node_count} nodes ({time.time() - t0:.0f} s)", flush=True)
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(rpm), n_modes=12)
    print(f"{key}: modes {np.round(freqs, 4).tolist()} ({time.time() - t0:.0f} s)", flush=True)
    n = mesh.node_count
    full = np.zeros((len(freqs), n, 6))
    full[:, free // 6, free % 6] = modes
    ring_id = np.full(n, -1)
    w = np.zeros(n)
    for k, (_z, nodes, wk) in enumerate(ring_profile(mesh, -1e9)):
        ring_id[nodes] = k
        w[nodes] = wk
    ids = np.empty(n, dtype=int)
    for nid, i in mesh.node_id_to_index.items():
        ids[i] = nid
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / f"shell-{key}.npz", coords=mesh.coords_array, node_ids=ids, freqs=freqs, modes=full,
             ring_id=ring_id, ring_w=w, box=box_node_mask(mesh), rpm=rpm)


def _fmt(a):
    return " ".join(f"{x:5.3f}" for x in a)


def step_compare():
    lines = ["# T8 RESULT: shell vs BeamDyn p20 mode shapes (es 0.25 mid skin; BeamDyn p20 from the existing .lin)", "",
             ("Generated by `tools/l3_mode_shapes.py compare`. Vectors: bending [u_x(z_i); u_y(z_i)], twist [theta_z(z_i)], "
             "combined [u_x; u_y; theta c/2] (c from the official yaml chord at z/117). Shell ring means arc-weighted "
             "(`ring_profile`), twist from `rigid_rotation_fit` on all ring nodes (A) or box nodes (B), linear in z "
             "between rings. Pairing: Hungarian on the combined MAC of shell modes 1..8 x BeamDyn 1..8 (square)."), ""]
    summary = {}
    for key in ("parked", "rated"):
        sh = np.load(OUT / f"shell-{key}.npz")
        bd = np.load(OUT / f"bd-{key}.npz")
        nb = N_MODES
        # BeamDyn: common phase from the combined complex vector, then real parts
        z, chord, (bux, buy, bth) = bd_station_modes(bd, nb)
        bd_f = bd["freqs"][:nb]
        bd_vec = {fit: vectors(bux, buy, bth, chord) for fit in ("all", "box")}
        n_sh = len(sh["freqs"])
        sh_f = sh["freqs"]
        res = {}
        for fit in ("all", "box"):
            sux, suy, sth, bad = shell_station_modes(sh["coords"], sh["modes"], sh["ring_id"], sh["ring_w"], sh["box"],
                                                     z, fit=fit)
            shv = vectors(sux, suy, sth, chord)
            shv = tuple(v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-300) for v in shv)
            bdv = tuple(v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-300) for v in bd_vec[fit])
            res[fit] = {"sh": (sux, suy, sth), "shv": shv, "bdv": bdv, "nan_modes": bad}
        labels_bd = family_labels(bux, buy, bth, chord)
        labels_sh = family_labels(*res["all"]["sh"], chord)
        lines += [f"## {key} ({'0' if key == 'parked' else RPM_RATED} rpm)", "",
                  "BeamDyn p20 modes (family by station energy, E = u_x, F = u_y, T = theta c/2): "
                  + ", ".join(f"{i + 1}:{labels_bd[i]} {bd_f[i]:.4f} Hz" for i in range(nb)), "",
                  "Shell modes: " + ", ".join(f"{i + 1}:{labels_sh[i]} {sh_f[i]:.4f} Hz" for i in range(n_sh)), ""]
        summary[key] = {"bd_freqs": bd_f.tolist(), "shell_freqs": sh_f.tolist(), "labels_bd": labels_bd,
                        "labels_sh": labels_sh}
        for fit, title in (("all", "A: twist from all ring nodes"), ("box", "B: twist from box nodes")):
            shv, bdv = res[fit]["shv"], res[fit]["bdv"]
            macs = [mac_matrix(s, b) for s, b in zip(shv, bdv, strict=True)]
            summary[key][fit] = {"mac_bend": macs[0].tolist(), "mac_twist": macs[1].tolist(),
                                 "mac_comb": macs[2].tolist(), "nan_modes": res[fit]["nan_modes"]}
            for tag, nsh in (("8x8", nb), ("12x8", n_sh)):
                pair, _ = pair_shell_to_bd(shv[2][:nsh], bdv[2])
                gap = (sh_f[pair] / bd_f - 1) * 100
                lines += [f"### {key}, {title}, pairing {tag} on the combined MAC", "",
                          "| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |",
                          "|---|---|---|---|---|---|---|---|"]
                for j in range(nb):
                    i = pair[j]
                    lines.append(f"| {j + 1} {labels_bd[j]} | {bd_f[j]:.4f} | {i + 1} {labels_sh[i]} | {sh_f[i]:.4f} | "
                                 f"{gap[j]:+.2f} | {macs[0][i, j]:.3f} | {macs[1][i, j]:.3f} | {macs[2][i, j]:.3f} |")
                summary[key][fit][f"pair_{tag}"] = pair.tolist()
                lines.append("")
            lines += [(f"Full MAC rows (shell modes 1..{n_sh} x BeamDyn 1..{nb}), {title}; shell modes within "
                       "3.7-4.9 Hz and all BeamDyn columns (bend | twist | comb):"), ""]
            for i in range(n_sh):
                if 3.7 <= sh_f[i] <= 4.9:
                    lines.append(f"- shell {i + 1} {labels_sh[i]} {sh_f[i]:.4f} Hz  bend [{_fmt(macs[0][i])}]  "
                                 f"twist [{_fmt(macs[1][i])}]  comb [{_fmt(macs[2][i])}]")
            lines += ["", "BeamDyn columns 1..8 are " + ", ".join(f"{j + 1}:{labels_bd[j]}" for j in range(nb)), ""]
            lines += ["BeamDyn modes within 3.7-4.9 Hz, column MAC over all shell modes (comb):", ""]
            for j in range(nb):
                if 3.7 <= bd_f[j] <= 4.9:
                    lines.append(f"- BD {j + 1} {labels_bd[j]} {bd_f[j]:.4f} Hz  comb over shell 1..{n_sh} [{_fmt(macs[2][:, j])}]")
            lines.append("")
    (OUT / "compare.json").write_text(json.dumps(summary, indent=1) + "\n")
    (OUT / "RESULT-tables.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", choices=["bd", "shell-parked", "shell-rated", "compare"])
    a = ap.parse_args()
    {"bd": step_bd, "shell-parked": lambda: step_shell("parked"), "shell-rated": lambda: step_shell("rated"),
     "compare": step_compare}[a.step]()


if __name__ == "__main__":
    main()
