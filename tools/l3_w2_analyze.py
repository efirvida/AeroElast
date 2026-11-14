#!/usr/bin/env python
"""L3 / T9 (W2): cap area, plateau thickness and end taper in a SONATA section dump (tools/l3_w2_sonata_cap.py npz).
args: npz matid t_yaml[m] hand_exact_per_cap[m2] ... prints per-cap area and loss in units of t^2 (round cutoff model: 0.8667 t^2 per cap)."""
import sys
import numpy as np
d = np.load(sys.argv[1], allow_pickle=True); mid = int(sys.argv[2]); t = float(sys.argv[3]); hand = [float(v) for v in sys.argv[4:]]
m = d["mat"] == mid; a = d["area"][m]; c = d["ctr"][m]
# separate the two caps by the sign of the centre y about the chord line is not valid when the airfoil is cambered: cluster in y gap
ys = np.sort(c[:, 1]); k = np.argmax(np.diff(ys)); thr = 0.5 * (ys[k] + ys[k + 1])
for name, sel in (("upper(y>thr)", c[:, 1] > thr), ("lower(y<thr)", c[:, 1] <= thr)):
    A = a[sel].sum(); x = c[sel, 0]
    print(f"{name}: n={sel.sum()} area={A:.5f} x[{x.min():.3f},{x.max():.3f}] extent={x.max()-x.min():.3f}")
    bins = np.linspace(x.min(), x.max() + 1e-9, 21); h = np.histogram(x, bins, weights=a[sel])[0] / np.diff(bins)
    print("   area density (m2/m) per x-bin (~vertical thickness):", np.round(h, 3).tolist())
    print(f"   plateau (median of central bins) {np.median(h[5:15]):.4f} m vs yaml t {t:.4f}")
if hand:
    A = a.sum(); H = 2 * np.mean(hand)
    print(f"total CarbonUD {A:.5f}; hand exact {sum(hand):.5f} -> loss {sum(hand)-A:.5f} m2 = {(sum(hand)-A)/t**2/2:.3f} t^2 per cap (round-cutoff model 0.867)")
