"""L3 T6 (W3) D2: does the station EI_edge gap explain the LC3 tip-load edgewise deflection excess at z 114.845 m?

Euler-Bernoulli cantilever, tip point load F: u(z*) = int_0^z* M(z) m(z) / EI(z) dz, M = F (L - z), m = z* - z (unit load at z*).
EI_deck(z) = deck K44 (edge, > K55 at every outboard station); EI_src(z) = r_src(eta) * EI_deck. Ratios r are EI_edge(source)/K44.
Shell thin-wall ratios come from logs/chord_scaled/w3/ei_edge.json; section-solver ratios are copied from the feature document
(L12 SONATA 2026; L16 SONATA 2021 and CalculiX solid). Outside the ratio table the ratio is held (variant H) or set to 1 (variant 1).
Check: u_deck vs the BeamDyn LC3 driver (31.904 m at node 56, gmsh-chord-scaled-mesh.md L65; 34.469 m at the tip).
usage: python tools/l3_w3_lc3_integral.py ei_edge.json
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from l3_w3_ei_edge import DECK, read_deck  # noqa: E402

js = json.load(open(sys.argv[1]))
L = js["L"]
rows = js["rows"]
eta_s = np.array([r["eta"] for r in rows])
src = {
    "shell mid (thin-wall)": (
        eta_s,
        np.array([r["mid"]["EI_edge"] / r["deck_edge"] for r in rows]),
    ),
    "shell OML (thin-wall)": (
        eta_s,
        np.array([r["oml"]["EI_edge"] / r["deck_edge"] for r in rows]),
    ),
    "shell mid A11* (reduced)": (
        eta_s,
        np.array([r["mid_red"]["EI_edge"] / r["deck_edge"] for r in rows]),
    ),
    "SONATA 2026 (L12)": (
        np.array([0.50, 0.75, 0.80, 0.85]),
        np.array([0.927, 0.861, 0.819, 0.746]),
    ),
    "SONATA 2021 (L16)": (
        np.array([0.30, 0.40, 0.50, 0.65, 0.75, 0.80, 0.85]),
        np.array([0.893, 0.945, 1.028, 0.910, 0.889, 0.807, 0.732]),
    ),
    "solid bonded (L16)": (np.array([0.50, 0.80, 0.85]), np.array([0.951, 0.839, 0.752])),
}
deck_eta, K = read_deck(DECK)
z = np.linspace(0, L, 117001)
eta = z / L
ei_d = np.interp(
    eta, deck_eta, K[:, 3, 3] - K[:, 2, 3] ** 2 / K[:, 2, 2]
)  # elastic-centre edge EI (K44 - K34^2/K33)
zs = 114.845
F = 1.0e6


def tip_u(ei, zstar):
    m = z <= zstar
    return np.trapezoid(F * (L - z[m]) * (zstar - z[m]) / ei[m], z[m])


ud = tip_u(ei_d, zs)
import yaml  # structural twist of the yaml (deg): the principal axes sit at this angle to the load direction x

_y = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "tests/IEA-15-240-RWT.yaml"))[
    "components"
]["blade"]["outer_shape_bem"]["twist"]
tw = np.radians(np.interp(eta, _y["grid"], np.degrees(_y["values"])))
ef = np.interp(eta, deck_eta, K[:, 4, 4] - K[:, 2, 4] ** 2 / K[:, 2, 2])  # deck flap K55


def tip_u_t(e_scale, f_scale, zstar):
    """as tip_u with twist coupling: x-compliance = cos^2/EI_edge + sin^2/EI_flap (principal axes rotated by the twist)."""
    c = np.cos(tw) ** 2 / (ei_d * e_scale) + np.sin(tw) ** 2 / (ef * f_scale)
    m = z <= zstar
    return np.trapezoid(F * (L - z[m]) * (zstar - z[m]) * c[m], z[m])


udt = tip_u_t(1.0, 1.0, zs)
print(
    "deck with twist coupling u(z*) = %.3f m (ratio to driver 31.904: %.3f)" % (udt, udt / 31.904)
)
print(
    "deck EB u(z*=114.845) = %.3f m (BeamDyn driver LC3: 31.904 m at node 56; ratio %.3f); tip %.3f m (driver 34.469)"
    % (ud, ud / 31.904, tip_u(ei_d, L) if False else np.trapezoid(F * (L - z) ** 2 / ei_d, z))
)
share = lambda a, b: (
    np.trapezoid(
        (F * (L - z) * (zs - z) / ei_d * (z <= zs))[(eta >= a) & (eta < b)],
        z[(eta >= a) & (eta < b)],
    )
    / ud
)  # noqa: E731
print(
    "share of u_deck by eta band: "
    + "  ".join(
        "%.2f-%.2f: %.1f%%" % (a, b, 100 * share(a, b))
        for a, b in ((0, 0.3), (0.3, 0.6), (0.6, 0.8), (0.8, 0.9), (0.9, 1.01))
    )
)
res = {}
for name, (e, r) in src.items():
    for var in ("H", "1"):
        lo, hi = (r[0], r[-1]) if var == "H" else (1.0, 1.0)
        rr = np.interp(eta, e, r, left=lo, right=hi)
        u = tip_u(ei_d / rr, zs) if False else tip_u(ei_d * rr, zs)
        res[(name, var)] = u / ud - 1
        print(
            "%-24s   with twist coupling (flap ratio 1): %+.1f %%"
            % ("", 100 * (tip_u_t(rr, 1.0, zs) / udt - 1))
        )
        print(
            "%-24s ratio table %.2f-%.2f  extension %s: u/u_deck - 1 = %+.1f %%   [shell FE LC3: +10.1 %%]"
            % (name, e[0], e[-1], var, 100 * (u / ud - 1))
        )
# the gap restricted to the shell thin-wall mid table only where it is measured (0.3-0.95)
json.dump(
    {f"{k[0]}|{k[1]}": v for k, v in res.items()},
    open(str(Path(sys.argv[1]).with_name("lc3_integral.json")), "w"),
    indent=1,
)
