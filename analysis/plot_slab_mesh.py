r"""What the layered cells actually look like, straight from the meshed deck.

Drawn from the .inp rather than from the geometry that was asked for, so what
is shown is what was solved: if the mesher lost a bridge or a layer came out a
different thickness, it appears here.

Four views:
 (a) the brine phase in 3D -- the spanning layers, with the ice bridges visible
     as holes through them.
 (b) a slice normal to the layers, showing the ice-brine-ice stack the load has
     to cross and how little of it is bridge.
 (c) a slice through one layer's plane, which is where the bridges live and the
     only place the transverse load path exists.
 (d) the element size distribution across the bridges, since the whole drained
     result rests on those being resolved.

    python3 analysis/plot_slab_mesh.py <Job-....inp> [out.png]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from audit_volume import parse_deck, expand


def load(path):
    nodes, elems, esets = parse_deck(path)
    keys = sorted(nodes)
    idx = {k: i for i, k in enumerate(keys)}
    P = np.asarray([nodes[k] for k in keys], dtype=float)
    name = next((n for n in ('Sphere_Only', 'SPHERE_ONLY') if n in esets), None)
    incl = [i for i in expand(esets[name]) if i in elems] if name else []
    conn_i = np.asarray([[idx[x] for x in elems[i]] for i in incl],
                        dtype=np.int64)
    return P, conn_i


def tet_faces(conn):
    f = np.vstack([conn[:, [0, 1, 2]], conn[:, [0, 1, 3]],
                   conn[:, [0, 2, 3]], conn[:, [1, 2, 3]]])
    f.sort(axis=1)
    _, i, c = np.unique(f, axis=0, return_index=True, return_counts=True)
    return f[i[c == 1]]          # surface = faces owned by one tet only


def panel_3d(ax, P, conn, L):
    surf = tet_faces(conn)
    # Draw the near side only: facets whose outward normal points toward the
    # viewer. Subsampling instead leaves a sparse dust of triangles that reads
    # as noise rather than as a surface, which is what an earlier version did.
    tris = P[surf]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    cen = tris.mean(axis=1)
    view = np.array([0.50, -0.72, 0.48])
    keep = np.einsum('ij,j->i', n, view) > 0
    tris = tris[keep]
    if len(tris) > 60000:
        tris = tris[:: max(1, len(tris) // 60000)]
    ax.add_collection3d(Poly3DCollection(
        tris, facecolor=fs.SKY, edgecolor='none', alpha=0.95))
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.view_init(elev=18, azim=-58)
    ax.set_title('(a) brine phase: spanning layers, bridges as holes')


def panel_slice(ax, P, conn, L, axis, at, title, xl, yl):
    """Brine elements whose centroid lies in a thin slab about `at`."""
    cen = P[conn].mean(axis=1)
    tol = 0.02 * L
    sel = np.abs(cen[:, axis] - at) < tol
    other = [i for i in (0, 1, 2) if i != axis]
    if sel.sum():
        ax.scatter(cen[sel][:, other[0]], cen[sel][:, other[1]],
                   s=5, color=fs.BLUE, alpha=0.55, linewidths=0)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_aspect('equal')
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title)


def panel_sizes(ax, P, conn, L):
    a, b, c, d = (P[conn[:, i]] for i in range(4))
    vol = np.abs(np.einsum('ij,ij->i', np.cross(b - a, c - a), d - a)) / 6.0
    h = (12.0 * vol / np.sqrt(2.0)) ** (1.0 / 3.0)     # edge of equivalent tet
    ax.hist(h / L, bins=45, color=fs.ORANGE, alpha=0.85)
    ax.axvline(np.median(h) / L, color=fs.BLACK, ls='--',
               label='median %.4f L' % (np.median(h) / L))
    ax.set_xlabel('element size / cell edge')
    ax.set_ylabel('brine elements')
    ax.set_title('(d) resolution of the brine phase')
    ax.legend(fontsize=11)


def main():
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'slab_mesh.png'
    P, conn = load(inp)
    if not len(conn):
        print('no inclusion elements in %s' % inp)
        return 1
    L = float(np.max(P))
    print('%s: %d nodes, %d brine tets, cell edge %.3f'
          % (os.path.basename(inp), len(P), len(conn), L))

    fig = plt.figure(figsize=(13.5, 11))
    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    panel_3d(ax1, P, conn, L)
    ax2 = fig.add_subplot(2, 2, 2)
    panel_slice(ax2, P, conn, L, 2, 0.5 * L,
                '(b) slice at z = L/2: the ice-brine stack', 'x', 'y')
    # a plane through the first layer, found from the brine centroids in x
    cen = P[conn].mean(axis=1)
    hist, edges = np.histogram(cen[:, 0], bins=60)
    xat = 0.5 * (edges[np.argmax(hist)] + edges[np.argmax(hist) + 1])
    ax3 = fig.add_subplot(2, 2, 3)
    panel_slice(ax3, P, conn, L, 0, xat,
                '(c) in the layer plane at x = %.3f: bridges are the gaps' % xat,
                'y', 'z')
    ax4 = fig.add_subplot(2, 2, 4)
    panel_sizes(ax4, P, conn, L)

    fig.tight_layout()
    fig.savefig(out, dpi=170)
    print('wrote %s' % out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
