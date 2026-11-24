"""Render the T4 meshes (uniform 0.25, N=10, N=20) at the root, mid span and tip as PNGs."""
import sys
from pathlib import Path

import pyvista as pv

pv.OFF_SCREEN = True
HERE = Path(__file__).parent
MESHES = [("uniform", "uniforme h 0.25"), ("n10", "N=10, AR 1"), ("n20", "N=20, AR 1")]
REGIONS = [("root", 0.0, 8.0), ("mid", 58.0, 62.0), ("tip", 107.0, 117.1), ("tipzoom", 113.0, 117.1)]


def region(mesh, z0, z1):
    z = mesh.cell_centers().points[:, 2]
    return mesh.extract_cells((z >= z0) & (z <= z1))


def main():
    grids = {k: pv.read(HERE / f"{k}.vtu") for k, _ in MESHES}
    for rname, z0, z1 in REGIONS:
        p = pv.Plotter(shape=(1, len(MESHES)), window_size=(1800, 700), off_screen=True)
        for j, (k, label) in enumerate(MESHES):
            p.subplot(0, j)
            sub = region(grids[k], z0, z1)
            skin = sub.threshold(0.5, scalars="is_web", invert=True)
            p.add_mesh(skin, color="lightgrey", show_edges=True, edge_color="black", line_width=0.6)
            p.add_text(f"{label}\nz {z0:g}-{z1:g} m, {grids[k].n_cells} elem", font_size=10)
            p.view_vector((0.0, 1.0, 0.0), viewup=(1.0, 0.0, 0.0))
            p.reset_camera()
        p.link_views()
        out = HERE / f"view_{rname}.png"
        p.screenshot(out)
        p.close()
        print(out)


if __name__ == "__main__":
    sys.exit(main())
