#!/usr/bin/env python3
"""Render the 3D matrix stress field of an RVE, for the critical-zone section.

The existing field render (viz/render_rve.py, mode `field`) draws the inclusions
as opaque grey bodies over the stress field, uses an absolute Pa scale, and
shows only the outer surface. All three work against the point the figure has to
make: the inclusions hide the matrix, the Pa scale means nothing without the
applied load, and the concentrations that matter are in the interior ligaments
between neighbouring channels, not on the boundary.

This renders the matrix alone. Inclusions are simply not drawn, so they read as
open voids and the field around them is visible; the stress is normalised by the
volume-weighted matrix mean, so the colour scale is a concentration factor and
carries no units; and the second panel cuts the cell open so the interior is
seen directly.

Two panels:
  (a) the intact cell, matrix only, viewed isometrically
  (b) the same cell cut on a mid-plane, exposing the ligaments

The colour scale is shared and clipped at a high percentile of the matrix field:
the extreme tail is a handful of tetrahedra and would otherwise flatten
everything else to one colour. The unclipped extreme is annotated instead.

    python3 render_stress_field.py base_z95.vtk stress_field_3d.png
"""
import sys

import numpy as np
import pyvista as pv

CLIP_PCT = 99.0          # colour-scale ceiling, as a percentile of the matrix field
# A map that is light at the bottom, not dark: half the matrix sits within a
# few percent of the mean, and on a dark-bottomed map like inferno that turns
# the whole cell black and hides the very concentrations the figure is for.
CMAP = 'YlOrRd'


def load_matrix(path):
    """Matrix cells only, with the field normalised to its volume-weighted mean."""
    mesh = pv.read(path)
    if 'material' not in mesh.cell_data:
        raise SystemExit('no material array: %s' % path)
    matrix = mesh.extract_cells(np.flatnonzero(mesh.cell_data['material'] == 0))
    mis = np.asarray(matrix.cell_data['mises'], float)
    vol = matrix.compute_cell_sizes(length=False, area=False,
                                    volume=True).cell_data['Volume']
    vol = np.abs(np.asarray(vol, float))
    mean = float((mis * vol).sum() / vol.sum())
    matrix.cell_data['scf'] = mis / mean
    return matrix, mean


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else 'base_z95.vtk'
    out = sys.argv[2] if len(sys.argv) > 2 else 'stress_field_3d.png'

    matrix, mean = load_matrix(src)
    scf = np.asarray(matrix.cell_data['scf'], float)
    hi = float(np.percentile(scf, CLIP_PCT))
    print('matrix cells %d   mean von Mises %.3e Pa' % (matrix.n_cells, mean))
    print('normalised field: p50 %.2f  p99 %.2f  max %.2f'
          % (np.percentile(scf, 50), hi, scf.max()))

    pv.global_theme.font.family = 'arial'
    pl = pv.Plotter(shape=(1, 2), window_size=(2200, 1150), off_screen=True,
                    border=False)
    sargs = dict(title='von Mises / matrix mean   ', vertical=False,
                 title_font_size=30, label_font_size=28, n_labels=5,
                 color='black', position_x=0.22, position_y=0.02,
                 width=0.56, height=0.055)
    outline = matrix.outline()
    b = matrix.bounds

    def frame(view_mesh, label, bar):
        # Render the bounding surface with backface culling. Drawing the volume
        # directly also draws the inward-facing skins of every void, which show
        # through the openings under different lighting and put a spurious grey
        # cast over parts of the cell; the field itself is uniform to within 1%
        # between the halves that appeared two-toned.
        surf = view_mesh.extract_surface()
        # Flat shading, no directional light. In an isometric view the top face
        # is lit and the two side faces sit in shadow, which on a light-based
        # colour map reads as a change in the field rather than in the lighting:
        # it put a grey cast over the side faces even though the field means of
        # the two halves agree to within 1%. Colour here therefore encodes the
        # field and nothing else, at the cost of the cube's solidity, which the
        # outline restores.
        pl.add_mesh(surf, scalars='scf', cmap=CMAP, clim=[1.0, hi],
                    show_scalar_bar=bar,
                    scalar_bar_args=sargs if bar else None,
                    lighting=False, backface_culling=True)
        pl.add_mesh(outline, color='#737373', line_width=2)
        pl.add_text(label, position=(0.02, 0.93), viewport=True,
                    font_size=22, color='black')
        pl.camera_position = 'iso'
        pl.camera.azimuth = 30
        pl.camera.elevation = 20
        pl.reset_camera()
        pl.camera.zoom(0.82)

    pl.subplot(0, 0)
    frame(matrix, '(a) intact cell', False)

    # A corner octant removed rather than a half: the cell still reads as a cube,
    # and the exposed faces cut through the interior ligaments.
    pl.subplot(0, 1)
    cx, cy = 0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])
    cut = (matrix.clip(normal='x', origin=(cx, 0, 0), invert=False)
                 .clip(normal='y', origin=(0, cy, 0), invert=False))
    keep = matrix.clip(normal='x', origin=(cx, 0, 0), invert=True)
    keep2 = matrix.clip(normal='x', origin=(cx, 0, 0), invert=False).clip(
        normal='y', origin=(0, cy, 0), invert=True)
    frame(keep + keep2, '(b) corner octant removed', True)

    pl.set_background('white')
    pl.screenshot(out, transparent_background=False)
    print('wrote %s  (scale clipped at %.2f, unclipped max %.1f)'
          % (out, hi, scf.max()))


if __name__ == '__main__':
    main()
