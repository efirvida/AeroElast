"""Close-up of the tip region of a mesh, both faces, coloured by element min angle."""
import sys
from pathlib import Path

import pyvista as pv

HERE = Path(__file__).parent
name = sys.argv[1] if len(sys.argv) > 1 else "n20"
z0, z1 = float(sys.argv[2]) if len(sys.argv) > 2 else 110.0, 117.1
m = pv.read(HERE / f"{name}.vtu")
z = m.cell_centers().points[:, 2]
sub = m.extract_cells((z >= z0) & (z <= z1))
sub = sub.threshold(0.5, scalars="is_web", invert=True)
p = pv.Plotter(shape=(2, 1), window_size=(2000, 1300), off_screen=True)
for i, v in enumerate([(0.0, 1.0, 0.0), (0.0, -1.0, 0.0)]):
    p.subplot(i, 0)
    p.add_mesh(sub, scalars="min_angle", cmap="RdYlGn", clim=[30, 90], show_edges=True, edge_color="black",
               line_width=0.5, scalar_bar_args={"title": "min angle (deg)"})
    p.add_text(f"{name} z {z0:g}-{z1:g} m, view {'+y' if i == 0 else '-y'}", font_size=11)
    p.view_vector(v, viewup=(1.0, 0.0, 0.0))
    p.reset_camera()
    p.camera.zoom(1.3)
out = HERE / f"close_{name}_{int(z0)}.png"
p.screenshot(out)
print(out)
