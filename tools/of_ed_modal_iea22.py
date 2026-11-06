"""ElastoDyn-only parked blade linearization for IEA-22 (analogue of the `ed_0` case of tools/of_bd_blade_modal.py).

python tools/of_ed_modal_iea22.py build
sbatch -p sequana_cpu_dev --ntasks=1 --time=15 --wrap "<binary> logs/check/of_ed/ed_0/case.fst"   # printed by build
python tools/of_ed_modal_iea22.py post

Inputs: the unmodified upstream IEA-22 ElastoDyn files pinned in tests/IEA22MW/reference/openfast (see the .source.json beside them).
The upstream .fst is the old AeroElasticSE format that OpenFAST v5 cannot read, so the v5-format primary file of the IEA-15 `ed_0`
case is used as the template for the simulation/linearization switches only (no turbine data lives there: every module but ElastoDyn
is off). Rotor parked (RotSpeed 0), gravity 0, only FlapDOF1/FlapDOF2/EdgeDOF on, as in the IEA-15 case. Outputs under logs/check/of_ed.
"""
import json, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import of_bd_blade_modal as m   # noqa: E402

REF = ROOT / 'tests/IEA22MW/reference/openfast'
INPUTS = ROOT / 'tests/IEA22MW/of_ed_modal'
LOGS = ROOT / 'logs/check/of_ed'


def build():
    case = INPUTS / 'ed_0'; out = LOGS / 'ed_0'
    case.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    fst = (ROOT / 'tests/IEA15MW/of_bd_modal/ed_0/case.fst').read_text()
    fst = m.replace(fst, 'EDFile', f'"{case / "ed.dat"}"')
    # The upstream ElastoDyn is the old AeroElasticSE layout (no PitchDOF, PtfmRefxt, PBrIner, ...) that v5 refuses. Convert by merging:
    # the IEA-15 `ed_0` file (v5 layout) is the template and every labelled parameter that exists exactly once in both files takes the
    # IEA-22 value. Parameters only in the template keep their IEA-15 value; none of them acts with all DOFs but the blade modes off.
    ed = m.empty_outputs((ROOT / 'tests/IEA15MW/of_bd_modal/ed_0/ed.dat').read_text())
    src = (REF / 'IEA-22-280-RWT_ElastoDyn.dat').read_text()
    label = re.compile(r'^[ \t]*(\S+)[ \t]+(\w[\w()]*)[ \t]+-', re.M)
    values = {k: v for v, k in label.findall(src)}
    merged, kept = [], []
    for v, key in label.findall(ed):
        if key in values and len(re.findall(rf'^[ \t]*\S+[ \t]+{re.escape(key)}[ \t]+-', src, re.M)) == 1:
            ed = m.replace(ed, key, values[key]); merged.append(key)
        else:
            kept.append(key)
    (LOGS / 'ed_0').joinpath('merge.json').write_text(json.dumps({'from_iea22': merged, 'kept_from_iea15_template': kept}, indent=1))
    print(len(merged), 'parameters from IEA-22;', len(kept), 'kept from the template:', kept)
    for key in re.findall(r'^\s*(?:True|False)\s+(\w*DOF\w*)\s', ed, re.M):
        ed = m.replace(ed, key, 'True' if key in ('FlapDOF1', 'FlapDOF2', 'EdgeDOF') else 'False')
    for key, value in {'RotSpeed': 0.0, 'SumPrint': 'False', 'NTwGages': 1, 'NBlGages': 3,   # gauge counts must match the template's gauge lists (output only)
                        'TwrFile': f'"{REF / "IEA-22-280-RWT_ElastoDyn_tower.dat"}"'}.items():
        ed = m.replace(ed, key, value)
    for i in range(1, 4):
        ed = m.replace(ed, f'BlPitch({i})', 0)
        ed = m.replace(ed, f'BldFile{i}', f'"{REF / "IEA-22-280-RWT_ElastoDyn_blade.dat"}"')
    (case / 'ed.dat').write_text(ed); (case / 'case.fst').write_text(fst); (out / 'case.fst').write_text(fst)
    print(f'{m.BINARY} {out / "case.fst"}')


if __name__ == '__main__':
    if sys.argv[1] == 'build':
        build()
    else:
        m.INPUTS, m.LOGS = INPUTS, LOGS
        m.post('ed_0', 1)
