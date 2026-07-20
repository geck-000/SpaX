"""Render publication-quality RVE figures from a SpaX mesh/result (pyvista,
offscreen). Two modes:

  micro  <Job-*.inp> <out.png>        microstructure: brine CHANNELS (blue) vs
                                      POCKETS (orange) in a faint ice cube.
  field  <Job-*.vtk> <out.png>        von Mises on the (warped) RVE, viridis.

Colours: Okabe-Ito colourblind-safe pair for the two inclusion classes; a
perceptually-uniform sequential map (viridis) for the stress magnitude -- no
rainbow. Usage: python render_rve.py <mode> <in> <out.png> [az] [el]
"""
import sys
import numpy as np
import pyvista as pv

BLUE, ORANGE = '#0072B2', '#E69F00'      # channels, pockets (Okabe-Ito, CVD-safe)
ICE = '#e8f2f8'                          # faint ice-cube tint


# ---------------------------------------------------------------- .inp parsing
def parse_inp(path):
    nodes, tets, mat, section, cur = {}, [], [], None, 0
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        if s.startswith('*'):
            low = s.lower()
            if low.startswith('*node') and 'output' not in low:
                section = 'node'
            elif low.startswith('*element') and 'output' not in low:
                section, cur = 'elem', (0 if 'matrix_only' in low else 1)
            else:
                section = None
            continue
        p = s.split(',')
        if section == 'node' and len(p) >= 4:
            nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
        elif section == 'elem' and len(p) >= 5:
            tets.append((int(p[1]), int(p[2]), int(p[3]), int(p[4])))
            mat.append(cur)
    return nodes, tets, np.array(mat, int)


def grid_from(nodes, tets, mat):
    labels = sorted(nodes)
    idx = {l: i for i, l in enumerate(labels)}
    pts = np.array([nodes[l] for l in labels], float)
    cells = np.empty((len(tets), 5), np.int64)
    cells[:, 0] = 4
    for k, t in enumerate(tets):
        cells[k, 1:] = [idx[t[0]], idx[t[1]], idx[t[2]], idx[t[3]]]
    g = pv.UnstructuredGrid(cells.ravel(),
                            np.full(len(tets), pv.CellType.TETRA, np.uint8), pts)
    g.cell_data['material'] = mat
    return g


# --------------------------------------------------------------- shared setup
def _plotter():
    pv.OFF_SCREEN = True
    p = pv.Plotter(off_screen=True, window_size=(1500, 1400), lighting='light kit')
    p.set_background('white')
    try:
        p.enable_anti_aliasing('ssaa')
    except Exception:
        pass
    return p


def _frame(p, bounds, az, el, cube=True):
    # faint tinted faces + a crisp dark wireframe so the cube reads clearly
    if cube:
        box = pv.Box(bounds=bounds)
        p.add_mesh(box, color=ICE, opacity=0.05, lighting=False)
        p.add_mesh(box.extract_all_edges(), color='#22303a', line_width=3.5)
    p.camera_position = 'iso'
    p.camera.azimuth = az
    p.camera.elevation = el
    p.enable_parallel_projection()
    p.reset_camera()               # fit the whole cube in frame
    p.camera.zoom(0.82)            # leave a margin for the legend/axes overlays
    p.add_axes(line_width=4, labels_off=False, color='#22303a',
               xlabel='x', ylabel='y', zlabel='z')


# ------------------------------------------------------------------ micro mode
def render_micro(inp, out, az, el):
    nodes, tets, mat = parse_inp(inp)
    g = grid_from(nodes, tets, mat)
    L = g.bounds[5] - g.bounds[4]
    incl = g.extract_cells(np.where(mat == 1)[0])
    # split inclusions into connected bodies; a body spanning most of the height
    # is a (percolating) channel, the compact ones are pockets.
    conn = incl.connectivity()
    rid = conn.cell_data['RegionId']
    kind = np.zeros(incl.n_cells, int)          # 0 pocket, 1 channel
    ptsz = conn.points[:, 2]
    for r in np.unique(rid):
        cids = np.where(rid == r)[0]
        sub = conn.extract_cells(cids)
        if (sub.points[:, 2].max() - sub.points[:, 2].min()) > 0.55 * L:
            kind[cids] = 1
    p = _plotter()
    for k, col, lab in ((1, BLUE, 'channels'), (0, ORANGE, 'pockets')):
        sel = np.where(kind == k)[0]
        if len(sel):
            surf = incl.extract_cells(sel).extract_surface()
            p.add_mesh(surf, color=col, smooth_shading=True, specular=0.25,
                       specular_power=15, ambient=0.30, diffuse=0.85, label=lab)
    _frame(p, g.bounds, az, el)
    p.add_legend(bcolor='white', border=True, size=(0.22, 0.13), loc='upper left',
                 face='circle', font_family='arial')
    p.screenshot(out, scale=2)
    print('micro %s: %d channels-cells, %d pocket-cells' %
          (out, int((kind == 1).sum()), int((kind == 0).sum())))


# ------------------------------------------------------------------ field mode
def render_field(vtk, out, az, el):
    g = pv.read(vtk)
    ctr = np.array(g.center)
    # von Mises lives in the ICE MATRIX (the soft brine is ~0 and would just read
    # as dark blobs). Show the matrix stress, cut open to reveal the interior; the
    # inclusions become light-grey context so the concentration around them shows.
    matrix = g.threshold(0.5, scalars='material', invert=True)   # material == 0
    incl = g.threshold(0.5, scalars='material')                  # material == 1
    matrix = matrix.cell_data_to_point_data()                    # smooth (not faceted)
    if 'U' in g.array_names:
        matrix = matrix.warp_by_vector('U', factor=3.0)
        incl = incl.warp_by_vector('U', factor=3.0)
    clip_n, clip_o = (0, 1, 0), ctr                              # cut on the y mid-plane
    mclip = matrix.clip(normal=clip_n, origin=clip_o)
    iclip = incl.clip(normal=clip_n, origin=clip_o)

    p = _plotter()
    p.add_mesh(iclip.extract_surface(), color='#c9ccd1', opacity=0.45,
               smooth_shading=True, specular=0.1, show_scalar_bar=False)
    p.add_mesh(mclip, scalars='mises', cmap='viridis', smooth_shading=True,
               specular=0.12, ambient=0.35, diffuse=0.85,
               scalar_bar_args=dict(title='von Mises  (Pa)', title_font_size=32,
                                    label_font_size=24, n_labels=5, fmt='%.1e',
                                    position_x=0.28, position_y=0.05, width=0.44,
                                    height=0.045, color='#22303a', font_family='arial'))
    _frame(p, g.bounds, az, el, cube=False)
    p.screenshot(out, scale=2)
    print('field %s written (matrix von Mises, y-cut)' % out)


def main():
    mode, src, out = sys.argv[1], sys.argv[2], sys.argv[3]
    az = float(sys.argv[4]) if len(sys.argv) > 4 else 40.0
    el = float(sys.argv[5]) if len(sys.argv) > 5 else 22.0
    (render_micro if mode == 'micro' else render_field)(src, out, az, el)


if __name__ == '__main__':
    main()
