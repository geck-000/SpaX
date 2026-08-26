r"""Figure 2: the three inclusion morphologies, as meshed.

Pockets, channels and layers are not alternatives; they are the succession the
material passes through with depth. Cold ice holds isolated pockets. Once the
brine percolates the channels appear. At the warm base the brine occupies
planes between ice platelets and the ice left between them carries the load
across those planes.

Every panel is read back from a solved deck rather than drawn from the geometry
that was requested, so what is shown is what was meshed: if the mesher lost an
inclusion or a layer came out a different thickness, it appears here.

Top row is the brine phase in three dimensions. Bottom row is a slice, chosen
per morphology to show the thing that matters -- for pockets and channels a
horizontal cut through the cell, for layers a cut in the plane of a layer,
where the ice bridges appear as gaps and the transverse load path is visible.

    python3 analysis/plot_fig2_morphologies.py <outdir>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from plot_slab_mesh import load, tet_faces

TMP = r'C:/Users/stirpeg2/.claude/jobs/06dae8ab/tmp'
CASES = [
    ('(a)', os.path.join(TMP, 'out_fig2', 'Job-FIG2pocket-utx.inp'),
     'cold ice: isolated inclusions', 2, None),
    ('(b)', os.path.join(TMP, 'out_fig2', 'Job-FIG2chan-utx.inp'),
     'once the brine percolates', 2, None),
    ('(c)',
     os.path.join(TMP, 'out_bracket', 'Job-BRK_p150_und-utx.inp'),
     'the warm base: brine in planes', 0, 'auto'),
]
VIEW = np.array([0.50, -0.72, 0.48])


def panel_3d(ax, P, conn, L, title, sub):
    surf = tet_faces(conn)
    tris = P[surf]
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    tris = tris[np.einsum('ij,j->i', n, VIEW) > 0]
    if len(tris) > 45000:
        tris = tris[:: max(1, len(tris) // 45000)]
    ax.add_collection3d(Poly3DCollection(tris, facecolor=fs.SKY,
                                         edgecolor='none', alpha=0.95))
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.view_init(elev=18, azim=-58)
    ax.text2D(0.02, 0.94, title, transform=ax.transAxes, fontsize=13,
              fontweight='bold', va='top')


def panel_slice(ax, P, conn, L, axis, at, xl, yl, note):
    cen = P[conn].mean(axis=1)
    sel = np.abs(cen[:, axis] - at) < 0.02 * L
    other = [i for i in (0, 1, 2) if i != axis]
    if sel.sum():
        ax.scatter(cen[sel][:, other[0]], cen[sel][:, other[1]], s=4,
                   color=fs.BLUE, alpha=0.55, linewidths=0)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(xl, fontsize=11); ax.set_ylabel(yl, fontsize=11)
    ax.text(0.03, 0.03, note, transform=ax.transAxes, fontsize=10.5,
            color='0.25')


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    fig = plt.figure(figsize=(13.2, 8.4))
    for i, (title, path, sub, axis, mode) in enumerate(CASES):
        if not os.path.exists(path):
            print('  missing %s' % path)
            continue
        P, conn = load(path)
        L = float(np.max(P))
        print('  %-24s %6d brine tets, L = %.3f' % (title, len(conn), L))

        ax = fig.add_subplot(2, 3, i + 1, projection='3d')
        panel_3d(ax, P, conn, L, title, sub)

        ax2 = fig.add_subplot(2, 3, i + 4)
        if mode == 'auto':
            # cut in the plane of a layer: the densest brine plane along x
            cen = P[conn].mean(axis=1)
            h, e = np.histogram(cen[:, 0], bins=60)
            at = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
            panel_slice(ax2, P, conn, L, 0, at, '$y$', '$z$',
                        'in the layer plane:\nbridges are the gaps')
        else:
            panel_slice(ax2, P, conn, L, axis, 0.5 * L, '$x$', '$y$',
                        'horizontal cut at $z=L/2$')

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = os.path.join(outdir, 'fig2_morphologies.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
