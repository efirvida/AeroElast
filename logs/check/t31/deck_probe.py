import sys, re, os, collections, numpy as np
tree, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, tree)
import tests.test_blade_ccx_parity as T
from aeroelast.core.mesh.io.writers import write_ccx_mesh
os.makedirs(out, exist_ok=True)
mesh, props, coords = T._blade()
load = T.CASES["tip_axial_tension"][0]
inp = os.path.join(out, "d.inp")
write_ccx_mesh(mesh, inp, properties=props, boundary_nodeset="RootNodes", solver_type="LinearStatic",
    load_nodeset="tip", load_vector=[float(x) for x in load]+[0,0,0], span_direction=(0,0,1.0), quadratic=True)
files = [inp] + [os.path.join(out, f) for f in os.listdir(out) if f.endswith(".msh") or f.endswith(".nam")]
txt = "".join(open(f).read() for f in sorted(set(files)))
print("TREE", tree, "mesh nodes", len(mesh.nodes), "mesh elems", len(mesh.elements), "sets", len(mesh.element_sets), "props", len(props))
# deck elements
blocks = re.split(r"(?im)^(?=\*)", txt)
deck_el = {}; elset_members = collections.defaultdict(set)
for b in blocks:
    m = re.match(r"\*ELEMENT, TYPE=(\w+)", b, re.I)
    if m:
        for ln in b.splitlines()[1:]:
            ln = ln.strip().rstrip(",")
            if ln and not ln.startswith("*"): deck_el.setdefault(m.group(1), set()).add(int(ln.split(",")[0]))
    m = re.match(r"\*ELSET, ELSET=(\w+)", b, re.I)
    if m:
        for ln in b.splitlines()[1:]:
            for t in ln.replace(" ","").split(","):
                if t.isdigit(): elset_members[m.group(1)].add(int(t))
ndeck = sum(len(v) for v in deck_el.values())
print("deck elements", ndeck, {k: len(v) for k, v in deck_el.items()})
# mesh id -> deck label: writer renumbers? compare by count per set
mesh_sets = {n: len(s.element_ids) if hasattr(s, "element_ids") else len(s.elements) for n, s in mesh.element_sets.items()}
# materials
mats = {}
for m in re.finditer(r"\*MATERIAL, NAME=(\S+)\n\*ELASTIC, TYPE=ENGINEERING CONSTANTS\n([^\n]*)\n([^\n]*)\n", txt):
    mats[m.group(1)] = ([float(x) for x in m.group(2).split(",")], None)
def eng(b): return b
secs = {}
for b in blocks:
    m = re.match(r"\*SHELL SECTION, ELSET=(\w+), COMPOSITE", b, re.I)
    if m:
        ps = []
        for ln in b.splitlines()[1:]:
            f = [x.strip() for x in ln.split(",")]
            if len(f) >= 3 and f[0]: ps.append((float(f[0]), f[2], f[3] if len(f) > 3 else ""))
        secs[m.group(1)] = ps
rows = []
iso=[n for n,p in props.items() if isinstance(p,dict)]
print('isotropic dict props:',len(iso),iso[:5])
for name, pr in props.items():
    if isinstance(pr,dict): continue
    el = "E" + name.upper(); ps = secs.get(el)
    a_ae = np.asarray(pr.a_matrix())[0, 0]; t_ae = pr.total_thickness
    if ps is None:
        rows.append((name, mesh_sets.get(name, 0), len(elset_members.get(el, ())), t_ae, np.nan, a_ae, np.nan)); continue
    t_d = sum(p[0] for p in ps); a_d = 0.0
    for t, mn, ori in ps:
        v, w = eng(mats[mn]); e1, e2, e3, n12, n13, n23, g12, g13 = v
        n21 = n12 * e2 / e1
        a_d += t * e1 / (1 - n12 * n21)   # angle 0 about span; Q11
    rows.append((name, mesh_sets.get(name, 0), len(elset_members.get(el, ())), t_ae, t_d, a_ae, a_d))
R = rows
miss = [r for r in R if np.isnan(r[4])]
print("sets with no composite section in deck:", len(miss), [(r[0], r[1]) for r in miss][:10])
cnt_mism = [r for r in R if r[1] != r[2]]
print("sets whose deck ELSET count != mesh set count:", len(cnt_mism), [(r[0], r[1], r[2]) for r in cnt_mism][:10])
ok = [r for r in R if not np.isnan(r[4])]
dt = np.array([r[4]/r[3]-1 for r in ok]); da = np.array([r[6]/r[5]-1 for r in ok])
print("thickness deck/AE-1: max|.| %.3e;  A11 deck/AE-1: median %.3e p95 %.3e max %.3e  (n=%d)" % (np.abs(dt).max(), np.median(da), np.percentile(np.abs(da),95), np.abs(da).max(), len(ok)))
# area-weighted total A11 (sum A11*area) both ways and per-span band
idx = mesh.node_id_to_index
def area(e):
    p = coords[[idx[n] for n in e.node_ids[:4]]]
    return 0.5*np.linalg.norm(np.cross(p[2]-p[0], p[3]-p[1])) if len(p)==4 else 0.5*np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0]))
byid = {e.id: e for e in mesh.elements}
tot_ae = tot_d = 0.0; band = collections.defaultdict(lambda: [0.0, 0.0])
pa = {r[0]: r for r in R}
for sn, s in mesh.element_sets.items():
    if sn not in pa or np.isnan(pa[sn][4]): continue
    for eid in (s.element_ids if hasattr(s, "element_ids") else [e.id for e in s.elements.values()]):
        e = byid[eid]; ar = area(e); z = coords[[idx[n] for n in e.node_ids]][:, 2].mean()
        b = int(z // 20) * 20
        band[b][0] += pa[sn][5]*ar; band[b][1] += pa[sn][6]*ar
        tot_ae += pa[sn][5]*ar; tot_d += pa[sn][6]*ar
print("sum(A11*area) deck/AE = %.6f" % (tot_d/tot_ae))
for b in sorted(band): print("  z %3d-%3d deck/AE %.6f  EA_ring(AE)=sum(A11*area)/20m = %.4e N" % (b, b+20, band[b][1]/band[b][0], band[b][0]/20.0))
# worst sets
for r in sorted(ok, key=lambda r: -abs(r[6]/r[5]-1))[:5]: print("  worst", r[0], "thk AE/deck %.6g/%.6g A11 AE/deck %.4e/%.4e" % (r[3], r[4], r[5], r[6]))
# BC and load
def card(name):
    for b in blocks:
        if b.upper().startswith(name): return b
for nm in ("*BOUNDARY", "*CLOAD"):
    b = card(nm); print(nm, "->", " | ".join(b.splitlines()[:4]) if b else None)
# nodes in root/tip nsets: mesh vs deck
rn = len(mesh.get_node_set("RootNodes").nodes); tn = len(mesh.get_node_set("tip").nodes)
def nset(n):
    for b in blocks:
        m = re.match(r"\*NSET, NSET=%s\b" % n, b, re.I)
        if m: return sum(1 for t in b.replace("\n", ",").split(",")[1:] if t.strip().isdigit())
print("root nodes mesh %d deck NROOTNODES %s ; tip nodes mesh %d deck NTIP %s" % (rn, nset("NROOTNODES"), tn, nset("NTIP")))
cl = card("*CLOAD")
if cl:
    ls = [l for l in cl.splitlines()[1:] if l.strip()]
    z = [l for l in ls if l.split(",")[1].strip() == "3"]
    print("CLOAD lines", len(ls), "dof3 lines", len(z), "sum fz", sum(float(l.split(",")[2]) for l in z))
zs = coords[:, 2]; print("z range mesh %.3f..%.3f" % (zs.min(), zs.max()), "tip z-extent", np.ptp([coords[idx[n.id], 2] for n in mesh.get_node_set("tip").nodes.values()]))
