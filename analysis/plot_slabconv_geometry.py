"""Render the slabconv mesh-convergence cell with pyvista.

A brine slab x in [0.4, 0.6] pierced by a single square ice bridge (y, z in
[0.4, 0.6]) -- the BRKB geometry in miniature, on the structured mesh that
resolves every phase boundary exactly at n % 10 == 0.  Ice is translucent so
the embedded brine slab (and the bridge hole through it) is visible.

    python3 analysis/plot_slabconv_geometry.py [out_png] [n]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'elements_ccx', 'tests'))
import smoothing_proto as S  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else 'slabconv_geometry.png'
N = int(sys.argv[2]) if len(sys.argv) > 2 else 24


def bridged(x, y, z, blo, bhi):
    if not (blo <= x < bhi):
        return 0
    return 0 if (0.4 <= y < 0.6) and (0.4 <= z < 0.6) else 1


nodes, tets, mat = S.mesh_box(N, 0.4, 0.6, 0.0, geom=bridged)

import pyvista as pv  # noqa: E402

grid = pv.UnstructuredGrid(
    {pv.CellType.TETRA: tets}, np.asarray(nodes, dtype=float))
grid.cell_data['mat'] = mat.astype(np.int32)

brine = grid.extract_cells(np.flatnonzero(mat == 1))
ice = grid.extract_cells(np.flatnonzero(mat == 0))

p = pv.Plotter(window_size=(1400, 1050), off_screen=True)
p.set_background('white')
# ice (matrix) translucent so the slab and bridge read through it
p.add_mesh(ice, color='#E8F2F8', opacity=0.25, show_edges=False,
           name='ice')
# brine (inclusion) solid
p.add_mesh(brine, color='#E69F00', opacity=1.0, show_edges=False,
           name='brine')
p.add_mesh(ice.extract_cells(ice['mat'] == 0), style='wireframe',
           color='#B0C4D8', line_width=0.5, opacity=0.35, name='wire')
p.camera_position = 'iso'
p.camera.azimuth = 35
p.camera.elevation = 25
p.add_axes(xlabel='x (load)', ylabel='y', zlabel='z')
p.save_graphic(OUT[:-4] + '.svg', title='slabconv cell')  # vector copy
p.screenshot(OUT)
p.close()
print('wrote %s and %s' % (OUT, OUT[:-4] + '.svg'))
