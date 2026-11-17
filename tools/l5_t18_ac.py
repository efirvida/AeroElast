"""T18: aerodynamic-centre fraction per station, production windIO blade_aero vs the AeroDyn-built one of the T3 tool."""
import json, sys
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "tools"))
import l5_oneway_bd_parity as par  # noqa: E402
import l5_t18_loads as lo  # noqa: E402
from aeroelast.models.blade.aerodynamics import load_blade_aero  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

cfg = lo.production_cfg(); b = cfg["bem"]
w = load_blade_aero(cfg["blade_file"], default_re=b["default_re"], neuralfoil_model=b["neuralfoil_model"], hub_radius=b["hub_radius"],
                    n_blades=b["n_blades"], viterna_ar=b["viterna_ar"], viterna_confidence_threshold=b["viterna_confidence_threshold"])
a = build_blade_aero_from_aerodyn(par.AD_PRIMARY)
def tab(ba): return {"r": [s.r for s in ba.stations], "ac": [s.airfoil.aerodynamic_center for s in ba.stations],
                    "pitch_axis": [s.pitch_axis for s in ba.stations], "chord": [s.chord for s in ba.stations], "twist_deg": [float(np.degrees(s.twist)) for s in ba.stations]}
out = {"windio": tab(w), "aerodyn": tab(a)}
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
for k, t in out.items():
    ac = np.array(t["ac"]); print(k, "n", len(ac), "ac min/max", ac.min(), ac.max(), "unique", np.unique(np.round(ac, 3))[:12])
# interpolate onto common r
rr = np.linspace(10, 118, 10)
for k, t in out.items(): print(k, np.round(np.interp(rr, t["r"], t["ac"]), 3), np.round(np.interp(rr, t["r"], t["chord"]), 2), np.round(np.interp(rr, t["r"], t["twist_deg"]), 2))
