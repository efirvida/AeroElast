"""T27: per-region A66 (mesher vs yaml) and thin/thick-wall Bredt GJ at one mid-span station.
argv: r/R (default 0.70). Needs mesh_head.pkl (mk_mesh.py), arcs.pkl (mk_arcs.py)."""
import sys, pickle, numpy as np, scipy.linalg as sl
sys.path.insert(0, "/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/gj")
import yamlref
from aeroelast.models.blade.model import build_rust_properties
W = "/scratch/leahk/eduardo.donestevez/fem-shell/logs/check/"
d = pickle.load(open(W + "gj/mesh_head.pkl", "rb")); ar = pickle.load(open(W + "gj/arcs.pkl", "rb"))
rr = float(sys.argv[1]) if len(sys.argv) > 1 else 0.70
cap = pickle.load(open(W + "input-fix/capture-T15m.pkl", "rb"))["fixtures"]
S1 = cap[[k for k in cap if k.startswith("s1_results")][0]]
zs = d["zs"]; ist = int(np.argmin(abs(zs / zs.max() - rr))); z = zs[ist]
L_blade = zs.max()
gj_fe = float(np.interp(z, S1["z"], S1["gj"])); gj_bd = float(np.interp(z, S1["z"], S1["gj_bd"]))
print(f"station idx {ist} z={z:.3f} r/R={z/L_blade:.4f}  FE S-1 GJ={gj_fe:.4e}  BD GJ={gj_bd:.4e}  FE/BD-1={gj_fe/gj_bd-1:+.2%}")

mats = {m["name"]: m for m in d["materials"]}
secs = {s["elementSet"]: s for s in d["numad_sections"]}
coords = d["coords"]; tang = d["tangs"][ist]; nodes = d["mem"][d["sts"][ist]]["nodes"]
p0 = coords[nodes].mean(axis=0); sd = (coords - p0) @ tang
a = np.array([1.0, 0, 0]); a = a - (a @ tang) * tang; a /= np.linalg.norm(a); b = np.cross(tang, a)
P2 = lambda x: np.array([(x - p0) @ a, (x - p0) @ b])
e2set = {}
for name, ids in d["esets"].items():
    if name in secs:
        for e in ids: e2set.setdefault(e, name)
web = set(d["esets"]["allShearWebEls"])

lam_cache = {}
def lam(name):
    if name not in lam_cache:
        p = build_rust_properties({"materials": d["materials"], "sections": [secs[name]]})[name]
        A = np.array(p.a_matrix()); D = np.array(p.d_matrix())
        lam_cache[name] = (A, D, float(p.total_thickness))
    return lam_cache[name]

verts = {}; segs = []   # seg: dict(ka,kb,name,isweb)
for eid, nl in d["elems"]:
    s = sd[nl]
    if not (s.min() < 0 <= s.max()): continue
    pts = []; n = len(nl)
    for k in range(n):
        i, j = nl[k], nl[(k + 1) % n]
        if (sd[i] < 0) != (sd[j] < 0):
            key = (min(i, j), max(i, j)); t = sd[i] / (sd[i] - sd[j])
            pts.append((key, P2(coords[i] + t * (coords[j] - coords[i]))))
    if len(pts) != 2: continue
    segs.append(dict(ka=pts[0][0], kb=pts[1][0], name=e2set[eid], isweb=eid in web))
    verts[pts[0][0]] = pts[0][1]; verts[pts[1][0]] = pts[1][1]
print("segments", len(segs), "vertices", len(verts), "distinct sections", len({s['name'] for s in segs}))

def seglen(s): return np.linalg.norm(verts[s["kb"]] - verts[s["ka"]])
def bredt(gt, V=None, extra=None):
    """Multi-cell Bredt-Batho: minimise complementary energy sum q^2 L/gt, Kirchhoff, torque 1. gt per segment (N/m)."""
    V = V or verts
    vk = list(V); vi = {k: i for i, k in enumerate(vk)}; ns = len(segs)
    Ls = np.zeros(ns); C = np.zeros(ns); M = np.zeros(ns); B = np.zeros((len(vk), ns))
    for i, s in enumerate(segs):
        pa, pb = V[s["ka"]], V[s["kb"]]
        Ls[i] = np.linalg.norm(pb - pa); C[i] = Ls[i] / gt[i]
        M[i] = pa[0] * pb[1] - pa[1] * pb[0]; B[vi[s["ka"]], i] -= 1; B[vi[s["kb"]], i] += 1
    N = sl.null_space(B); m = N.T @ M; Cm = N.T @ (C[:, None] * N)
    y = np.linalg.solve(Cm, m); q = N @ y; q /= (M @ q)
    return 1.0 / np.sum(q * q * C)

def a66(A): return A[2, 2], 1.0 / np.linalg.inv(A)[2, 2]
gt_ours = np.array([a66(lam(s["name"])[0])[1] for s in segs])
GJ_ours = bredt(gt_ours)

# ---------------- yaml reference ----------------
eta = z / L_blade
YL = yamlref.layers(eta)
cap_w = [(L["ex"]["width"], L["ex"]["end_nd_arc"] - L["ex"]["start_nd_arc"]) for L in YL if L["name"] == "Spar_Cap_SS"][0]
P = cap_w[0] / cap_w[1]
for L in YL:
    ex = L["ex"]; s0, s1 = ex.get("start_nd_arc"), ex.get("end_nd_arc")
    if L["name"].startswith("TE_reinforcement_SS"): s0, s1 = 0.0, ex["width"] / P
    if L["name"].startswith("TE_reinforcement_PS"): s1, s0 = 1.0, 1.0 - ex["width"] / P
    L["s0"], L["s1"] = s0, s1
sLE = [0.5 * (L["ex"]["start_nd_arc"] + L["ex"]["end_nd_arc"]) for L in YL if L["name"] == "LE_reinforcement"][0]
print(f"yaml: eta={eta:.4f} perimeter P={P:.3f} m (from cap width), s_LE={sLE:.4f}")

def ref_plies(s):
    out = []
    for L in YL:
        if L["web"] is None and L["s0"] is not None and L["s0"] - 1e-9 <= s <= L["s1"] + 1e-9 and L["t"] > 0:
            out.append([L["mat"], L["t"], L["ang"]])
    return out
def ref_web(wname):
    return [[L["mat"], L["t"], L["ang"]] for L in YL if L["web"] == wname and L["t"] > 0]

def ply_Q(mat):
    el = mat["elastic"]; E = el["E"]
    if isinstance(E, (list, tuple)):
        E1, E2 = float(E[0]), float(E[1]); nu12 = float(el["nu"][0]); G12 = float(el["G"][0]) if el.get("G") else E1 / (2 * (1 + nu12))
    else:
        E1 = E2 = float(E); nu12 = float(el["nu"]); G12 = float(el["G"]) if el.get("G") else E1 / (2 * (1 + nu12))
    nu21 = nu12 * E2 / E1; den = 1 - nu12 * nu21
    return np.array([[E1 / den, nu12 * E2 / den, 0], [nu12 * E2 / den, E2 / den, 0], [0, 0, G12]])
def q66_ply(m, ang):
    Q = ply_Q(mats[m]); th = np.deg2rad(ang); c, s_ = np.cos(th), np.sin(th)
    T = np.array([[c*c, s_*s_, 2*c*s_], [s_*s_, c*c, -2*c*s_], [-c*s_, c*s_, c*c - s_*s_]]); R = np.diag([1, 1, 2.0])
    return (np.linalg.inv(T) @ Q @ R @ T @ np.linalg.inv(R))[2, 2]
def lam_of(plies):
    p = build_rust_properties({"materials": d["materials"], "sections": [dict(type="shell", layup=plies, elementSet="x")]})["x"]
    return np.array(p.a_matrix()), np.array(p.d_matrix()), float(p.total_thickness)

# region arcs (m from LE; HP<0). Linear map arc->nd_arc per side, fitted on the cap edges.
ka = ar["key_arcs"][:, int(np.argmin(abs(ar["ispan"] - z)))]
ca = {L["name"]: L for L in YL}
def fit(capname, a_edges):
    L = ca[capname]; s_e = np.array([L["s0"], L["s1"]]); return np.polyfit(a_edges, s_e, 1)
fit_ss = fit("Spar_Cap_SS", [ka[8], ka[9]][::-1] if False else [ka[8], ka[9]]) if False else None
# mesher caps: LP keys 8..9 = 'b','c' (0.645,1.545), HP keys 3..4 = 'c','b' (-1.516,-0.616)
fl_ss = np.polyfit([ka[8], ka[9]], [ca["Spar_Cap_SS"]["s1"], ca["Spar_Cap_SS"]["s0"]], 1)
fl_ps = np.polyfit([ka[3], ka[4]], [ca["Spar_Cap_PS"]["s1"], ca["Spar_Cap_PS"]["s0"]], 1)
to_s = lambda am: np.polyval(fl_ss if am > 0 else fl_ps, am)
regs = {}   # j -> name
for s in segs:
    if not s["isweb"]:
        j = int(s["name"][:2]); regs[j] = s["name"].split("_", 2)[2].rsplit("_", 1)[0] if s["name"].endswith(("_A", "_B", "_C")) else s["name"].split("_", 2)[2]
def region_ref(j, nsamp=21):
    lo, hi = sorted([ka[j], ka[j + 1]]); 
    if hi - lo < 1e-6: return None
    As = []; hs = []
    for am in np.linspace(lo, hi, nsamp + 2)[1:-1]:
        pl = ref_plies(to_s(am))
        if not pl: continue
        A, Dm, h = lam_of(pl); As.append(a66(A)); hs.append(h)
    if not As: return None
    return np.mean(As, axis=0), np.mean(hs), np.min([x[1] for x in As]), np.max([x[1] for x in As])

# ---------------- per-region table ----------------
rows = []
for j in sorted(regs):
    names = sorted({s["name"] for s in segs if not s["isweb"] and int(s["name"][:2]) == j})
    for nm in names:
        Ao, Do, ho = lam(nm); ln = sum(seglen(s) for s in segs if s["name"] == nm)
        rows.append((j, nm, ln, ho, *a66(Ao)))
print("\nPer-region in-plane shear stiffness A66 (N/m).  A66 = A[2,2]; Aeff = 1/(A^-1)66 (Bredt-relevant). ref = yaml ply stack")
print(f"{'sect':22s} {'j':>2s} {'Lseg m':>7s} | {'h ours mm':>9s} {'A66 ours':>10s} {'Aeff ours':>10s} | {'h yaml mm':>9s} {'A66 yaml':>10s} {'Aeff yaml':>10s} [min,max over region] | {'A66 d%':>7s} {'Aeff d%':>8s}")
gt_ref = []; tab = {}
refcache = {}
for j, nm, ln, ho, Ao66, Aoe in rows:
    r = region_ref(j)
    if r is None: continue
    refcache[j] = r
    (Ar66, Are), hr, mn, mx = r
    tab[nm] = (Ar66, Are)
    print(f"{nm:22s} {j:2d} {ln:7.3f} | {ho*1e3:9.2f} {Ao66:10.4e} {Aoe:10.4e} | {hr*1e3:9.2f} {Ar66:10.4e} {Are:10.4e} [{mn:.3e},{mx:.3e}] | {Ao66/Ar66-1:+7.1%} {Aoe/Are-1:+8.1%}")
# webs
wnames = sorted({s["name"] for s in segs if s["isweb"]})
wref = {}
for k, nm in enumerate(wnames):
    Ao, Do, ho = lam(nm); ln = sum(seglen(s) for s in segs if s["name"] == nm)
    cand = "web0" if nm.lower().find("0") >= 0 and nm.lower().find("1") < 0 else None
    wref[nm] = (ho, ln, a66(Ao))
print("\nweb sections at the station:", {n: (round(v[0] * 1e3, 2), round(v[1], 3)) for n, v in wref.items()})
for wn in ("web0", "web1"):
    A, Dm, h = lam_of(ref_web(wn)); print(f"yaml {wn}: h={h*1e3:.2f} mm A66={a66(A)[0]:.4e} Aeff={a66(A)[1]:.4e}")
pickle.dump(dict(wref=wref, tab=tab), open(W + "gj/tab_%03d.pkl" % round(rr * 100), "wb"))

# ---------------- Bredt with reference laminate ----------------
def ref_gt(s):
    if s["isweb"]:
        wn = "web0" if verts[s["ka"]][0] < 0 else "web1"      # placeholder; replaced below by side
        return None
    j = int(s["name"][:2]); r = refcache.get(j)
    return r[0][1] if r else None
# web side: closest web ordering by chordwise coordinate
webcx = {}
for nm in wnames:
    xs = [0.5 * (verts[s["ka"]][0] + verts[s["kb"]][0]) for s in segs if s["name"] == nm]
    webcx[nm] = np.mean(xs)
order = sorted(webcx, key=webcx.get)
print("web order by x:", order)
# web0 offset_y_pa -0.35 (LE side?) -- compare thicknesses instead of guessing
gt_ref = []
for s, g in zip(segs, gt_ours):
    if s["isweb"]:
        ho = lam(s["name"])[2]
        cands = {wn: lam_of(ref_web(wn)) for wn in ("web0", "web1")}
        wn = min(cands, key=lambda w: abs(cands[w][2] - ho))
        gt_ref.append(a66(cands[wn][0])[1])
    else:
        r = refcache.get(int(s["name"][:2])); gt_ref.append(r[0][1] if r else g)
gt_ref = np.array(gt_ref)
GJ_ref = bredt(gt_ref)
print(f"\nThin-wall Bredt (multi-cell, midline geometry of the mid-surface mesh):\n  ours (mesher laminate): {GJ_ours:.4e}  ({GJ_ours/gj_bd-1:+.2%} vs BD)\n  yaml laminate A66:      {GJ_ref:.4e}  ({GJ_ref/gj_bd-1:+.2%} vs BD)   ref/ours-1={GJ_ref/GJ_ours-1:+.2%}")
print(f"  FE S-1 {gj_fe:.4e} ({gj_fe/gj_bd-1:+.2%} vs BD); FE/Bredt-ours-1={gj_fe/GJ_ours-1:+.2%}")

# ---------------- thick-wall variants ----------------
V = {k: np.array(v, float) for k, v in verts.items()}
cen = np.mean(list(V.values()), axis=0)
nb = {k: [] for k in V}
for i, s in enumerate(segs):
    if not s["isweb"]: nb[s["ka"]].append(i); nb[s["kb"]].append(i)
def seg_normal(i):
    s = segs[i]; tv = V[s["kb"]] - V[s["ka"]]; n = np.array([tv[1], -tv[0]]); n /= np.linalg.norm(n)
    mid = 0.5 * (V[s["ka"]] + V[s["kb"]])
    return n if n @ (mid - cen) > 0 else -n      # outward (skins only; approx. star-shaped about centroid)
def shift_variant(fn):
    """fn(seg_index) -> signed shift along outward normal (+ = toward OML). Vertex = mean over adjacent skin segs."""
    off = {}
    for k, lst in nb.items():
        if not lst: continue
        off[k] = np.mean([seg_normal(i) * fn(i) for i in lst], axis=0)
    return {k: V[k] + off.get(k, 0) for k in V}
def tot_t(i): return lam(segs[i]["name"])[2]
def gel_sign(i):
    lay = secs[segs[i]["name"]]["layup"]; mats_ = [x[0] for x in lay]
    if "Gelcoat" in mats_ and mats_[0] == "Gelcoat": return +1   # list order OML -> IML
    if "Gelcoat" in mats_ and mats_[-1] == "Gelcoat": return -1  # list order IML -> OML
    return +1
def centroid_shift(i):
    """Through-thickness centroid of the ply shear stiffness (Q66-weighted), signed toward the OML."""
    nm = segs[i]["name"]; lay = secs[nm]["layup"]; sg = gel_sign(i)
    h = np.array([x[1] for x in lay]); zc = np.cumsum(h) - 0.5 * h - 0.5 * h.sum()   # list-order coordinate, 0 = midplane
    w = []
    for (m, t, ang) in lay:
        w.append(q66_ply(m, ang) * t)
    w = np.array(w); zbar = float((w * zc).sum() / w.sum())
    return sg * (-zbar)      # list coordinate grows toward the IML when the list runs OML->IML
GJ = {}
GJ["mid"] = bredt(gt_ours)
GJ["OML"] = bredt(gt_ours, shift_variant(lambda i: +0.5 * tot_t(i)))
GJ["IML"] = bredt(gt_ours, shift_variant(lambda i: -0.5 * tot_t(i)))
GJ["Q66-centroid"] = bredt(gt_ours, shift_variant(centroid_shift))
shifts = np.array([centroid_shift(i) for i in range(len(segs)) if not segs[i]["isweb"]])
Lsk = np.array([seglen(segs[i]) for i in range(len(segs)) if not segs[i]["isweb"]])
print("\nThick-wall variants (laminate = mesher's):")
for k, v in GJ.items(): print(f"  closed-section Bredt, wall path = {k:13s}: {v:.4e}  vs mid {v/GJ['mid']-1:+.2%}  vs BD {v/gj_bd-1:+.2%}")
print(f"  Q66-centroid shift: length-weighted mean {np.sum(shifts*Lsk)/Lsk.sum()*1e3:+.2f} mm (+ = toward OML), range [{shifts.min()*1e3:+.1f},{shifts.max()*1e3:+.1f}] mm; mean skin h={np.mean([tot_t(i) for i in range(len(segs)) if not segs[i]['isweb']])*1e3:.1f} mm")
# open-section St Venant term: sum 4 D66 L over every wall (skin + web) -- the shell's D66 twisting stiffness the thin-wall Bredt omits
def D66_eff(nm):
    A, D, h = lam(nm); return D[2, 2]
open_term = sum(4.0 * D66_eff(s["name"]) * seglen(s) for s in segs)
open_skin = sum(4.0 * D66_eff(s["name"]) * seglen(s) for s in segs if not s["isweb"])
print(f"  open-section term sum 4*D66*L: all walls {open_term:.4e} ({open_term/GJ['mid']:+.2%} of Bredt mid; skins only {open_skin/GJ['mid']:+.2%})")
GJ_thick = GJ["mid"] + open_term
print(f"  thick-wall GJ = Bredt(mid) + open term = {GJ_thick:.4e}  vs BD {GJ_thick/gj_bd-1:+.2%}")

# what would close the gap: uniform wall-path shift toward the OML (fraction f of each segment's half thickness) and uniform G scale
from scipy.optimize import brentq
try: f_req = brentq(lambda f: bredt(gt_ours, shift_variant(lambda i: f * 0.5 * tot_t(i))) - gj_bd, -1.5, 2.5)
except ValueError: f_req = float('nan')
print(f"  closing the Bredt gap to BD needs the wall path at {f_req:+.2f} x (h/2) toward the OML = {f_req*0.5*np.mean([tot_t(i) for i in range(len(segs)) if not segs[i]['isweb']])*1e3:+.1f} mm at the mean skin h; or all-wall G x{gj_bd/GJ['mid']:.3f}")
# points of BD
pb = lambda v: 100 * (v - gj_bd) / gj_bd
print(f"  [pts of BD] FE {pb(gj_fe):+.2f} | Bredt thin ours {pb(GJ_ours):+.2f} | Bredt thin yaml-laminate {pb(GJ_ref):+.2f} | thick (mid+open) {pb(GJ_thick):+.2f} | Q66-centroid path {pb(GJ['Q66-centroid']):+.2f}")
print(f"  [contributions, pts of BD] A66 laminate mismatch {pb(GJ_ref)-pb(GJ_ours):+.2f} | thick-wall open term {pb(GJ_thick)-pb(GJ_ours):+.2f} | wall-path centroid {pb(GJ['Q66-centroid'])-pb(GJ_ours):+.2f} | FE below its thin Bredt {pb(gj_fe)-pb(GJ_ours):+.2f} | unexplained after thick {-pb(GJ_thick):+.2f}")
# isotropic-circular sanity of the second-order wall-thickness correction: (1+(t/2Rm)^2)
Rm = np.sqrt(0.5 * abs(sum(0.5 * (V[s['ka']][0] * V[s['kb']][1] - V[s['ka']][1] * V[s['kb']][0]) * 2 for s in segs if not s['isweb'])) / np.pi) if False else None
res = dict(rr=rr, z=z, gj_fe=gj_fe, gj_bd=gj_bd, GJ_ours=GJ_ours, GJ_ref=GJ_ref, GJ=GJ, open_term=open_term, open_skin=open_skin)
pickle.dump(res, open(W + "gj/res_%03d.pkl" % round(rr * 100), "wb"))
