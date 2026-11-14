"""Span-wise survey of a mesh: one PNG per z segment, both skin faces and the webs, coloured by min angle."""
import sys
from pathlib import Path

import pyvista as pv

HERE = Path(__file__).parent
# Fractions of the blade length, so every blade is surveyed alike (IEA-15: 0-10, 10-22, ... 112-117 m).
FRACTIONS = [0.0, 0.085, 0.19, 0.31, 0.445, 0.6, 0.75, 0.87, 0.957, 1.0001]


def main(name):
    m = pv.read(HERE / f"{name}.vtu")
    zc = m.cell_centers().points[:, 2]
    out_dir = HERE / f"survey_{name}"
    out_dir.mkdir(exist_ok=True)
    zlo, zhi = m.bounds[4], m.bounds[5]
    edges = [round(zlo + f * (zhi - zlo), 1) for f in FRACTIONS]
    for z0, z1 in zip(edges[:-1], edges[1:]):
        sub = m.extract_cells((zc >= z0) & (zc <= z1))
        skin = sub.threshold(0.5, scalars="is_web", invert=True)
        webs = sub.threshold(0.5, scalars="is_web")
        p = pv.Plotter(shape=(3, 1), window_size=(2200, 1650), off_screen=True)
        views = [((0.0, 1.0, 0.0), skin, "+y skin"), ((0.0, -1.0, 0.0), skin, "-y skin"), ((-1.0, 0.0, 0.0), webs, "webs (-x)")]
        for i, (v, part, label) in enumerate(views):
            p.subplot(i, 0)
            if part.n_cells:
                p.add_mesh(part, scalars="min_angle", cmap="RdYlGn", clim=[30, 90], show_edges=True,
                           edge_color="black", line_width=0.4, show_scalar_bar=(i == 0),
                           scalar_bar_args={"title": "min angle (deg)"})
            p.add_text(f"{name} z {z0:g}-{z1:g} m  {label}", font_size=11)
            p.view_vector(v, viewup=(1.0, 0.0, 0.0) if i < 2 else (0.0, 1.0, 0.0))
            p.reset_camera()
            p.camera.zoom(1.25)
        out = out_dir / f"z{z0:05.1f}-{z1:05.1f}.png"
        p.screenshot(out)
        p.close()
        print(out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "n20")
