#!/usr/bin/env python3
"""Re-point the legacy ``BladeMesh`` / ``RotorMesh`` case configs at the windIO 2.x official inputs.

gmsh-official-mesher T6b (S6, user decision L-S41).  After T6a/T6c ``BladeMesh`` and ``RotorMesh`` build with gmsh and
reject a windIO 1.x yaml, an ``excel_file`` and the NuMAD-only keys, so every legacy config fails to load.  This tool edits
ONLY the mesh generator block (``mesh.generator.type`` and ``mesh.generator.params``) of each config, line by line, so that
comments, key order and every other key stay byte-identical (PyYAML is the only YAML library of the environment and it
drops comments, so the file is never re-dumped).  The result is re-parsed and compared with the planned generator block;
a mismatch skips the file instead of writing it.

Rules (one per mapping):

* ``type: BladeMesh`` -> ``type: GmshBladeMesh``; ``type: RotorMesh`` keeps its type.
* ``yaml_file`` naming a windIO 1.x file (no ``windIO_version``, or 1.x) of a known turbine -> the relative path, from the
  config's directory, to that turbine's upstream v2 file.  A 1.x file of no known turbine (for example the rigid-material
  variant of the IEA-15) has no v2 twin: the config is skipped and listed.
* the NuMAD-only keys ``n_samples``, ``airfoil_spacing``, ``span_grading``, ``airfoil_dir``, ``refine_tip``,
  ``max_bay_length`` and ``laminate_per_row`` are dropped.
* ``element_size`` within 0.15-0.25 m is kept; absent, or above 0.25 m (AGENTS.md rule 9), it is replaced by the reference
  mesh ``chord_divisions: 20`` / ``span_aspect: 1`` (added only when the config carries no ``chord_divisions``).  Values
  below 0.15 m are kept and reported.
* ``skin_reference`` is written explicitly for ``BladeMesh`` configs with the value the runner resolved
  (``_blade_skin_reference``: the key, else ``"mid"``), so the skin choice never rests on a default.  ``RotorMesh`` has no such
  key and is left alone.
* ``blade_file`` (top level, under ``aero:``, any nesting) naming a windIO 1.x yaml with a v2 twin is re-pointed the same way
  (S8); a config whose generator or ``blade_file`` reads a 1.x yaml WITHOUT a twin is ``git mv``-ed to
  ``tests/archive/windio1_configs/<path under tests>`` unless something else names it by file name.
* a config with ``excel_file`` (the UTD xlsx blade) is not rewritten: it is ``git mv``-ed to
  ``tests/archive/numad_xlsx_configs/<path under tests>``, with a README, unless something else names it by file name.

Other keys that still point at a windIO 1.x file (aero/BEM inputs...) are never touched, only reported.

Default is ``--dry-run`` (a per-file plan); ``--write`` applies it.  Idempotent: a migrated config yields an empty plan.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]

NUMAD_KEYS = (
    "n_samples", "airfoil_spacing", "span_grading", "airfoil_dir", "refine_tip", "max_bay_length", "laminate_per_row",
)
# file name of a windIO 1.x turbine -> its upstream v2 file (relative to the repository root)
UPSTREAM_V2 = {
    "IEA-15-240-RWT.yaml": "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml",
    "IEA-3.4-130-RWT.yaml": "tests/IEA34/reference/upstream/IEA-3.4-130-RWT.yaml",
    "IEA-22-280-RWT.yaml": "tests/IEA22MW/reference/upstream/IEA-22-280-RWT.yaml",
}
IEA15_V2 = UPSTREAM_V2["IEA-15-240-RWT.yaml"]
ES_MIN, ES_MAX = 0.15, 0.25
ARCHIVE = "tests/archive/numad_xlsx_configs"
ARCHIVE_W1 = "tests/archive/windio1_configs"
_BLADE_FILE_RE = re.compile(r"^(\s*)blade_file:\s*([\"']?)([^\"'#\s]+)\2(.*)$")
REFERENCE_CHORD_DIVISIONS, REFERENCE_SPAN_ASPECT = 20, 1
_TEXT_SUFFIXES = {".py", ".sh", ".srm", ".yaml", ".yml", ".md", ".xml", ".toml", ".json", ".txt", ".sbatch", ".slurm", ".cfg"}
_GEN_RE = re.compile(r"type:\s*[\"']?(BladeMesh|RotorMesh)\b")
_README_W1 = """# windio1_configs

Case configs that read a windIO 1.x blade yaml with no upstream v2 twin (the rigid-material copy
`tests/IEA-15-240-RWT_rigid.yaml`), archived by gmsh-official-mesher T6b (S8, user decision L-S46 of
`odd/tasks/gmsh-official-mesher.md`: windIO 1.x is retired, the work continues on v2 only).  `BladeMesh` / `RotorMesh` reject a 1.x
yaml, so these configs no longer load.  Their relative paths (`../../...`) are relative to their former directory under `tests/`.
"""
_ARCHIVE_README = """# numad_xlsx_configs

Case configs that built the IEA-15 blade from the UTD NuMAD spreadsheet (`excel_file`), archived by
gmsh-official-mesher T6b (S6, user decision L-S41 of `odd/tasks/gmsh-official-mesher.md`).

The NuMAD mesher is retired (S5) and `BladeMesh` / `RotorMesh` reject `excel_file` (S7a): these configs no longer load.  The
UTD blade has a thinner outboard trailing-edge reinforcement than the official windIO input (AGENTS.md rule 7,
`docs/blade_input_divergence_utd_vs_official.md`), so none of them is a current anchor.  They are kept as a record of the
input-sensitivity runs only; their relative paths (`../../...`) are relative to their former directory under `tests/`.
"""


# ----------------------------------------------------------------------------------------------------------------------
# plan
# ----------------------------------------------------------------------------------------------------------------------
@dataclass
class Plan:
    path: Path
    action: str = "unchanged"  # rewritten | archived | skipped | unchanged
    changes: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)  # items to list in the report (element_size < 0.15, other turbines...)
    reason: str = ""
    new_text: str | None = None
    target: Path | None = None  # archive destination
    readme: str | None = None  # README text for the archive directory of the target


def windio_major(path: Path) -> int | None:
    """Major windIO version of a file; 1 when it declares none; None when unreadable."""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    m = re.search(r"^windIO_version:\s*['\"]?(\d+)", text, re.MULTILINE)
    return int(m.group(1)) if m else 1


def _indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _is_content(line: str) -> bool:
    s = line.strip()
    return bool(s) and not s.startswith("#")


def _block_end(lines: list[str], start: int, indent: int) -> int:
    """Index after the block that starts below ``lines[start]`` (a line with ``indent``): first content line at indent <= it."""
    i = start + 1
    while i < len(lines):
        if _is_content(lines[i]) and _indent(lines[i]) <= indent:
            break
        i += 1
    # trailing blanks/comments belong to what follows, not to the block
    while i > start + 1 and not _is_content(lines[i - 1]):
        i -= 1
    return i


def _locate(lines: list[str]) -> tuple[int, int, int, int, int]:
    """(type_line, params_line, params_end, child_indent, generator_indent) of ``mesh.generator``; ValueError if not found."""
    mesh = next((i for i, l in enumerate(lines) if re.match(r"mesh:\s*(#.*)?$", l)), None)
    if mesh is None:
        raise ValueError("no top-level 'mesh:' block")
    mesh_end = _block_end(lines, mesh, 0)
    gen = next((i for i in range(mesh + 1, mesh_end) if re.match(r"\s+generator:\s*(#.*)?$", lines[i])), None)
    if gen is None:
        raise ValueError("no 'generator:' block under mesh")
    gen_end = _block_end(lines, gen, _indent(lines[gen]))
    typ = next((i for i in range(gen + 1, gen_end) if re.match(r"\s+type:\s*\S", lines[i])), None)
    par = next((i for i in range(gen + 1, gen_end) if re.match(r"\s+params:\s*(#.*)?$", lines[i])), None)
    if typ is None or par is None:
        raise ValueError("generator block without a block-style 'type:' and 'params:'")
    par_end = _block_end(lines, par, _indent(lines[par]))
    first = next((i for i in range(par + 1, par_end) if _is_content(lines[i])), None)
    if first is None:
        raise ValueError("empty params block")
    return typ, par, par_end, _indent(lines[first]), _indent(lines[gen])


def _units(lines: list[str], par: int, par_end: int, child: int) -> list[tuple[str, int, int]]:
    """(key, first_line, end_line) of every ``key:`` entry of the params block, with its deeper continuation lines."""
    out: list[tuple[str, int, int]] = []
    i = par + 1
    while i < par_end:
        m = re.match(r"([A-Za-z0-9_]+)\s*:", lines[i][child:]) if _is_content(lines[i]) and _indent(lines[i]) == child else None
        if not m:
            i += 1
            continue
        j = i + 1
        while j < par_end and _is_content(lines[j]) and _indent(lines[j]) > child:
            j += 1
        out.append((m.group(1), i, j))
        i = j
    return out


def _rel(target: Path, start: Path) -> str:
    return Path(os.path.relpath(target, start)).as_posix()


def _other_v1_references(data: Any, cfg: Path, root: Path) -> list[str]:
    """Scalar values outside ``mesh.generator`` that name an existing windIO 1.x turbine file (reported, never edited)."""
    found: list[str] = []

    def walk(node: Any, trail: str) -> None:
        if trail == "mesh.generator":
            return
        if isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{trail}.{k}" if trail else str(k))
        elif isinstance(node, list):
            for n, v in enumerate(node):
                walk(v, f"{trail}[{n}]")
        elif isinstance(node, str) and node.endswith((".yaml", ".yml")) and "RWT" in Path(node).name and not trail.endswith("blade_file"):
            target = (cfg.parent / node).resolve()
            if target.is_file() and windio_major(target) == 1:
                found.append(f"{trail}: {node}")

    walk(data, "")
    return found


def _plan_generator(cfg: Path, root: Path, text: str) -> Plan | None:
    """The plan for the ``mesh.generator`` block, or None when it is not a BladeMesh/RotorMesh one."""
    plan = Plan(cfg)
    try:
        data = yaml.safe_load(text)
        gen = data["mesh"]["generator"]
        gtype, params = gen["type"], gen.get("params") or {}
    except Exception as exc:  # unparsable or not the generator layout
        plan.action, plan.reason = "skipped", f"cannot read mesh.generator ({exc.__class__.__name__})"
        return plan
    if gtype not in ("BladeMesh", "RotorMesh"):
        return None
    plan.notes += [f"other v1 reference (not edited) {r}" for r in _other_v1_references(data, cfg, root)]

    if params.get("excel_file"):
        return _plan_archive(plan, cfg, root)

    yaml_file = params.get("yaml_file")
    if not yaml_file:
        plan.action, plan.reason = "skipped", "no yaml_file in the generator params"
        return plan
    src = (cfg.parent / yaml_file).resolve()
    major = windio_major(src)
    new_yaml: str | None = None
    if major is None:
        plan.action, plan.reason = "skipped", f"yaml_file {yaml_file} not found"
        return plan
    if major == 1:
        v2 = UPSTREAM_V2.get(src.name)
        if v2 is None:
            return _plan_archive(plan, cfg, root, ARCHIVE_W1, f"windIO 1.x yaml {yaml_file} has no upstream v2 twin", _README_W1)
        new_yaml = _rel((root / v2).resolve(), cfg.parent.resolve())
        plan.changes.append(f"yaml_file: {yaml_file} -> {new_yaml}")
        if v2 != IEA15_V2:
            plan.notes.append(f"non-IEA-15 turbine: {src.name} -> {v2}")
    elif src.name != Path(IEA15_V2).name:
        plan.notes.append(f"non-IEA-15 turbine (already v2): {src.name}")

    try:
        lines = text.splitlines(keepends=True)
        typ, par, par_end, child, _ = _locate(lines)
    except ValueError as exc:
        plan.action, plan.reason = "skipped", str(exc)
        return plan
    units = _units(lines, par, par_end, child)
    present = {k for k, _, _ in units}
    pad = " " * child
    edits: dict[int, str | None] = {}  # line index -> replacement (None deletes); continuation lines are deleted too
    inserts: list[tuple[int, str]] = []  # (before line index, text)

    expect_type = gtype
    if gtype == "BladeMesh":
        expect_type = "GmshBladeMesh"
        edits[typ] = re.sub(r"(type:\s*[\"']?)BladeMesh\b", r"\1GmshBladeMesh", lines[typ], count=1)
        plan.changes.append("type: BladeMesh -> GmshBladeMesh")

    expect = dict(params)
    if new_yaml is not None:
        i = next(i for k, i, _ in units if k == "yaml_file")
        edits[i] = re.sub(r"^(\s*yaml_file:\s*)([\"']?)[^\"'#\s]+\2", lambda m: f"{m.group(1)}{m.group(2)}{new_yaml}{m.group(2)}", lines[i], count=1)
        expect["yaml_file"] = new_yaml
    for key, first, end in units:
        if key in NUMAD_KEYS:
            for n in range(first, end):
                edits[n] = None
            expect.pop(key)
            plan.changes.append(f"drop {key}")

    after_yaml = next(e for k, _, e in units if k == "yaml_file") if "yaml_file" in present else par + 1
    new_keys: list[str] = []
    if gtype == "BladeMesh" and "skin_reference" not in present:
        inserts.append((after_yaml, f"{pad}skin_reference: mid  # explicit: what the runner resolved for BladeMesh configs\n"))
        expect["skin_reference"] = "mid"
        new_keys.append("skin_reference: mid")

    es = params.get("element_size", None)
    sizing_at = after_yaml
    drop_es = False
    if "element_size" in present:
        try:
            es_value = float(es)
        except (TypeError, ValueError):
            plan.action, plan.reason = "skipped", f"element_size {es!r} is not a number"
            return plan
        if es_value < ES_MIN - 1e-9:
            plan.notes.append(f"element_size {es_value:g} < {ES_MIN} kept")
        elif es_value > ES_MAX + 1e-9:
            drop_es = True
    if drop_es or "element_size" not in present:
        if "element_size" in present:
            first, end = next((f, e) for k, f, e in units if k == "element_size")
            for n in range(first, end):
                edits[n] = None
            sizing_at = first
            expect.pop("element_size")
            plan.changes.append(f"drop element_size {es:g}")
        if "chord_divisions" not in present:
            new_chord = [f"{pad}chord_divisions: {REFERENCE_CHORD_DIVISIONS}  # reference mesh N20, span AR 1 (AGENTS.md rule 9)\n",
                         f"{pad}span_aspect: {REFERENCE_SPAN_ASPECT}\n"]
            inserts.append((sizing_at, "".join(new_chord)))
            expect["chord_divisions"], expect["span_aspect"] = REFERENCE_CHORD_DIVISIONS, REFERENCE_SPAN_ASPECT
            new_keys.append(f"chord_divisions: {REFERENCE_CHORD_DIVISIONS}, span_aspect: {REFERENCE_SPAN_ASPECT}")
    if new_keys:
        plan.changes.append("add " + "; ".join(new_keys))

    if not plan.changes:
        return plan
    out: list[str] = []
    for n, line in enumerate(lines):
        # inserts placed at the same index keep the order they were queued in
        out += [t for at, t in inserts if at == n]
        if n in edits:
            if edits[n] is not None:
                out.append(edits[n])  # type: ignore[arg-type]
        else:
            out.append(line)
    out += [t for at, t in inserts if at >= len(lines)]
    new_text = "".join(out)

    try:
        back = yaml.safe_load(new_text)["mesh"]["generator"]
        ok = back["type"] == expect_type and back["params"] == expect
        ok = ok and _without_generator(yaml.safe_load(new_text)) == _without_generator(data)
    except Exception:
        ok = False
    if not ok:
        plan.action, plan.reason, plan.changes = "skipped", "line edit does not reproduce the planned generator block", []
        return plan
    plan.action, plan.new_text = "rewritten", new_text
    return plan


def _blade_file_pass(cfg: Path, root: Path, text: str, plan: Plan) -> Plan:
    """Re-point every ``blade_file`` (any nesting) naming a windIO 1.x yaml with a v2 twin; archive when there is no twin."""
    lines = text.splitlines(keepends=True)
    out = list(lines)
    mapping: dict[str, str] = {}
    for n, line in enumerate(lines):
        m = _BLADE_FILE_RE.match(line.rstrip("\n"))
        if not m or not m.group(3).endswith((".yaml", ".yml")):
            continue
        value = m.group(3)
        src = (cfg.parent / value).resolve()
        if windio_major(src) != 1:
            continue
        v2 = UPSTREAM_V2.get(src.name)
        if v2 is None:
            return _plan_archive(plan, cfg, root, ARCHIVE_W1, f"blade_file {value} is windIO 1.x with no upstream v2 twin", _README_W1)
        new = _rel((root / v2).resolve(), cfg.parent.resolve())
        out[n] = f"{m.group(1)}blade_file: {m.group(2)}{new}{m.group(2)}{m.group(4)}" + ("\n" if line.endswith("\n") else "")
        mapping[value] = new
        plan.changes.append(f"blade_file: {value} -> {new}")
        if v2 != IEA15_V2:
            plan.notes.append(f"non-IEA-15 turbine: {src.name} -> {v2}")
    if not mapping:
        return plan
    new_text = "".join(out)

    def swap(node: Any) -> Any:
        if isinstance(node, dict):
            return {k: (mapping.get(v, v) if k == "blade_file" and isinstance(v, str) else swap(v)) for k, v in node.items()}
        if isinstance(node, list):
            return [swap(v) for v in node]
        return node

    try:
        ok = yaml.safe_load(new_text) == swap(yaml.safe_load(text))
    except Exception:
        ok = False
    if not ok:
        plan.action, plan.reason, plan.changes = "skipped", "line edit does not reproduce the planned blade_file change", []
        return plan
    plan.action, plan.new_text = "rewritten", new_text
    return plan


def plan_file(cfg: Path, root: Path) -> Plan | None:
    """The plan for one config (generator block, then ``blade_file``), or None when it has neither to migrate."""
    text = cfg.read_text(encoding="utf-8")
    has_gen = bool(_GEN_RE.search(text))
    if not has_gen and not any(_BLADE_FILE_RE.match(l.rstrip("\n")) for l in text.splitlines()):
        return None
    plan = _plan_generator(cfg, root, text) if has_gen else None
    if plan is not None and plan.action in ("skipped", "archived"):
        return plan
    base = plan or Plan(cfg)
    cur = base.new_text if base.new_text is not None else text
    result = _blade_file_pass(cfg, root, cur, base)
    if plan is None and result.action == "unchanged":
        return None
    return result


def _without_generator(data: dict) -> dict:
    out = dict(data)
    mesh = dict(out.get("mesh") or {})
    gen = dict(mesh.get("generator") or {})
    gen.pop("type", None)
    gen.pop("params", None)
    mesh["generator"] = gen
    out["mesh"] = mesh
    return out


def _referenced_elsewhere(cfg: Path, root: Path) -> list[str]:
    """Files under tests/tools/scripts/cases/src that name the config by file name (other than the config itself)."""
    name = cfg.name
    hits: list[str] = []
    for top in ("tests", "tools", "scripts", "cases", "src"):
        base = root / top
        if not base.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if d not in {".git", "__pycache__", "node_modules", "numad_xlsx_configs"}]
            for fn in filenames:
                p = Path(dirpath) / fn
                if p == cfg or p.suffix not in _TEXT_SUFFIXES or p.stat().st_size > 2_000_000:
                    continue
                try:
                    if name in p.read_text(encoding="utf-8", errors="ignore"):
                        hits.append(p.relative_to(root).as_posix())
                except OSError:
                    pass
    return hits


def _plan_archive(plan: Plan, cfg: Path, root: Path, base: str = ARCHIVE, why: str = "excel_file config",
                  readme: str = _ARCHIVE_README) -> Plan:
    rel = cfg.resolve().relative_to((root / "tests").resolve()) if cfg.resolve().is_relative_to((root / "tests").resolve()) else Path(cfg.name)
    refs = _referenced_elsewhere(cfg, root)
    if refs:
        plan.action = "skipped"
        plan.reason = f"{why} kept in place, named by: " + ", ".join(refs[:5])
        return plan
    plan.action = "archived"
    plan.target = root / base / rel
    plan.readme = readme
    plan.changes.append(f"{why}: git mv -> {plan.target.relative_to(root).as_posix()}")
    return plan


# ----------------------------------------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------------------------------------
def collect(paths: list[Path], root: Path) -> list[Plan]:
    files: list[Path] = []
    for p in paths:
        if p.is_file():
            files.append(p)
        else:
            files += [f for f in p.rglob("*") if f.suffix in (".yaml", ".yml") and f.is_file()]
    archive = (root / "tests/archive").resolve()
    plans = []
    for f in sorted(set(files)):
        if f.resolve().is_relative_to(archive):
            continue
        plan = plan_file(f, root)
        if plan is not None:
            plans.append(plan)
    return plans


def apply(plan: Plan, root: Path) -> None:
    if plan.action == "rewritten":
        assert plan.new_text is not None
        plan.path.write_text(plan.new_text, encoding="utf-8")
    elif plan.action == "archived":
        assert plan.target is not None
        plan.target.parent.mkdir(parents=True, exist_ok=True)
        moved = subprocess.run(["git", "-C", str(root), "mv", str(plan.path), str(plan.target)], capture_output=True, text=True)
        if moved.returncode != 0:  # not a git tree (or untracked): plain move
            shutil.move(str(plan.path), str(plan.target))
        base = next(root / b for b in (ARCHIVE, ARCHIVE_W1) if plan.target.is_relative_to(root / b))
        if not (base / "README.md").exists():
            (base / "README.md").write_text(plan.readme or _ARCHIVE_README, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("paths", nargs="*", type=Path, help="config files or directories (default: <root>/tests)")
    ap.add_argument("--root", type=Path, default=ROOT, help="repository root holding the upstream v2 files (default: this repo)")
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="print the plan only (default)")
    mode.add_argument("--write", action="store_true", help="apply the plan")
    args = ap.parse_args(argv)
    root = args.root.resolve()
    paths = [p.resolve() for p in (args.paths or [root / "tests"])]
    plans = collect(paths, root)

    counts: dict[str, int] = {}
    for pl in plans:
        counts[pl.action] = counts.get(pl.action, 0) + 1
        if pl.action == "unchanged" and not pl.notes:
            continue
        rel = pl.path.relative_to(root).as_posix() if pl.path.is_relative_to(root) else str(pl.path)
        print(f"{pl.action.upper():9s} {rel}" + (f"  [{pl.reason}]" if pl.reason else ""))
        for c in pl.changes:
            print(f"            {c}")
        for n in pl.notes:
            print(f"            note: {n}")
    print("summary: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())) + f" (total {len(plans)})"
          + ("  [dry-run: nothing written]" if not args.write else ""))
    if args.write:
        for pl in plans:
            apply(pl, root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
