"""Render publication-quality RVE figures from a SpaX mesh/result (pyvista,
offscreen). Three modes:

  micro  <Job-*.inp> <out.png>        microstructure: brine CHANNELS (blue) vs
                                      POCKETS (orange) in a faint ice cube.
  mesh   <Job-*.inp> <out.png>        the periodic FE mesh itself, one corner
                                      octant removed so the conforming
                                      matrix/inclusion interface shows.
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
MATRIX = '#dde6ee'                       # ice matrix in the cut-away mesh view
EDGE = '#7d8b99'                         # element edges


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


def _trim(path, pad=24):
    """Crop the uniform white margin left by the offscreen camera, keeping a
    small pad, so the figure carries no dead space into the paper."""
    try:
        from PIL import Image, ImageChops
    except ImportError:
        return
    im = Image.open(path).convert('RGB')
    # tolerance crop: the renderer leaves a scattering of near-white pixels that
    # would otherwise defeat an exact-white bounding box
    diff = ImageChops.difference(im, Image.new('RGB', im.size, (255, 255, 255)))
    bbox = diff.convert('L').point(lambda v: 255 if v > 8 else 0).getbbox()
    if not bbox:
        return
    x0, y0, x1, y1 = bbox
    im.crop((max(0, x0 - pad), max(0, y0 - pad),
             min(im.width, x1 + pad), min(im.height, y1 + pad))).save(path)


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
    if '--noaxes' not in sys.argv:
        p.add_axes(line_width=4, labels_off=False, color='#22303a',
                   xlabel='x', ylabel='y', zlabel='z')


def _classify(incl, L):
    """Split the inclusion cells into connected bodies; a body spanning most of
    the cell height is a (percolating) channel, the compact ones are pockets.
    Returns a per-cell array, 0 = pocket, 1 = channel."""
    conn = incl.connectivity()
    rid = conn.cell_data['RegionId']
    kind = np.zeros(incl.n_cells, int)
    for r in np.unique(rid):
        cids = np.where(rid == r)[0]
        sub = conn.extract_cells(cids)
        if (sub.points[:, 2].max() - sub.points[:, 2].min()) > 0.55 * L:
            kind[cids] = 1
    return kind


# ------------------------------------------------------------------ micro mode
def render_micro(inp, out, az, el):
    nodes, tets, mat = parse_inp(inp)
    g = grid_from(nodes, tets, mat)
    L = g.bounds[5] - g.bounds[4]
    incl = g.extract_cells(np.where(mat == 1)[0])
    kind = _classify(incl, L)
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


# ------------------------------------------------------------------- mesh mode
def _mesh_render(inp, out, az, el, cut):
    """The periodic tetrahedral mesh of a real generated deck. With cut=True one
    corner octant of cells is removed -- whole tets are dropped (cell-centre
    test) rather than sliced, so every element drawn is a genuine element of the
    solved model -- exposing the interior discretisation and the conforming
    matrix/inclusion interface. With cut=False the intact cell shows the strictly
    periodic surface mesh and the inclusion traces on the cell faces."""
    nodes, tets, mat = parse_inp(inp)
    g = grid_from(nodes, tets, mat)
    b = g.bounds
    L = b[5] - b[4]
    ctr = [(b[0] + b[1]) / 2, (b[2] + b[3]) / 2, (b[4] + b[5]) / 2]
    if cut:
        cc = g.cell_centers().points
        # octant facing the camera, so the cut opens toward the reader
        oct_ = (cc[:, 0] < ctr[0]) & (cc[:, 1] > ctr[1]) & (cc[:, 2] > ctr[2])
        keep = g.extract_cells(np.where(~oct_)[0])
    else:
        keep = g
    kmat = keep.cell_data['material']

    # near-flat shading: the cell faces must read as one material, not as three
    # brightnesses, so the mesh lines carry the 3D form instead of the lighting
    shade = dict(show_edges=True, edge_color=EDGE, line_width=0.35,
                 specular=0.0, ambient=0.72, diffuse=0.30)
    p = _plotter()
    matrix = keep.extract_cells(np.where(kmat == 0)[0]).extract_surface()
    p.add_mesh(matrix, color=MATRIX, label='ice matrix', **shade)
    incl = keep.extract_cells(np.where(kmat == 1)[0])
    if incl.n_cells:
        kind = _classify(incl, L)
        for k, col, lab in ((1, BLUE, 'brine channels'), (0, ORANGE, 'brine pockets')):
            sel = np.where(kind == k)[0]
            if len(sel):
                p.add_mesh(incl.extract_cells(sel).extract_surface(), color=col,
                           label=lab, **shade)
    _frame(p, g.bounds, az, el)
    if '--nolegend' not in sys.argv:      # the paper labels the phases in TikZ
        p.add_legend(bcolor='white', border=True, size=(0.235, 0.115),
                     loc='upper left', face='circle', font_family='arial')
    p.screenshot(out, scale=2)
    _trim(out)
    print('mesh %s: %d of %d tets shown (%d inclusion)' %
          (out, keep.n_cells, g.n_cells, int((kmat == 1).sum())))


def render_mesh(inp, out, az, el):
    _mesh_render(inp, out, az, el, cut=False)


def render_meshcut(inp, out, az, el):
    _mesh_render(inp, out, az, el, cut=True)


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
    {'micro': render_micro, 'mesh': render_mesh, 'meshcut': render_meshcut,
     'field': render_field}[mode](src, out, az, el)


if __name__ == '__main__':
    main()
