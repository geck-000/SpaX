r"""A single clean view of a layered cell, for use as a call-out.

plot_slab_mesh.py draws four diagnostic panels and takes the cell edge on
trust; this takes the edge from the node coordinates and draws one thing --
the brine phase in 3D, with the ice bridges appearing as holes through the
spanning layers. That is the whole point of the morphology and the only thing
a call-out needs to show.

    python3 analysis/plot_slab_cell.py <Job-....inp> [out.png]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

BRINE = '#E8A33D'


def read_inp(path):
    """Node coordinates and the element sets, straight from the deck."""
    nodes, elems, cur, setname = {}, {}, None, None
    with open(path, encoding='utf8', errors='replace') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('*'):
                low = s.lower()
                if low.startswith('*node'):
                    cur = 'node'
                elif low.startswith('*element'):
                    cur = 'elem'
                    setname = None
                    for part in s.split(','):
                        if 'elset' in part.lower():
                            setname = part.split('=')[1].strip()
                    elems.setdefault(setname, [])
                else:
                    cur = None
                continue
            if not s or cur is None:
                continue
            bits = [b.strip() for b in s.split(',') if b.strip()]
            if cur == 'node' and len(bits) >= 4:
                nodes[int(bits[0])] = (float(bits[1]), float(bits[2]),
                                       float(bits[3]))
            elif cur == 'elem' and len(bits) >= 5:
                elems[setname].append([int(b) for b in bits[1:5]])
    return nodes, elems


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    inp = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'slab_cell.png'

    nodes, elems = read_inp(inp)
    used = np.unique([n for k in elems for e in elems[k] for n in e])
    xyz = np.array([nodes[i] for i in used if i in nodes])
    L = float(xyz.max())                     # meshed nodes only, never the RPs

    brine_sets = [k for k in elems
                  if k and any(t in k.lower() for t in ('brine', 'incl',
                                                        'sphere', 'soft'))]
    if not brine_sets:                       # fall back to the smaller set
        brine_sets = [min(elems, key=lambda k: len(elems[k]))]
    idx = np.unique([n for k in brine_sets for e in elems[k] for n in e])
    pts = np.array([nodes[i] for i in idx])

    # Solid surfaces, not a point cloud: take the boundary of the brine
    # element set -- the faces that belong to exactly one tet -- so the layers
    # render as sheets and the bridges as holes through them.
    from collections import Counter
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    faces = Counter()
    for k in brine_sets:
        for e in elems[k]:
            for f in ((e[0], e[1], e[2]), (e[0], e[1], e[3]),
                      (e[0], e[2], e[3]), (e[1], e[2], e[3])):
                faces[tuple(sorted(f))] += 1
    hull = [f for f, c in faces.items() if c == 1]

    tris = np.array([[nodes[i] for i in f] for f in hull])
    # drop the faces lying on the cell boundary: they are cut surfaces, not
    # the brine-ice interface, and they hide the interior when shaded
    cen = tris.mean(axis=1)
    tol = 1e-6
    on_face = ((np.abs(cen) < tol) | (np.abs(cen - L) < tol)).any(axis=1)
    tris = tris[~on_face]

    fig = plt.figure(figsize=(4.6, 4.4))
    ax = fig.add_subplot(111, projection='3d')
    pc = Poly3DCollection(tris, facecolors=BRINE, edgecolors=BRINE,
                          linewidths=0.05, alpha=0.95)
    ax.add_collection3d(pc)

    for a, b in ((0, 0), (0, L), (L, 0), (L, L)):
        ax.plot([a, a], [b, b], [0, L], color='0.35', lw=0.7)
    for z in (0, L):
        ax.plot([0, L, L, 0, 0], [0, 0, L, L, 0], [z] * 5, color='0.35', lw=0.7)

    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-58)
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(out, dpi=200, transparent=True)
    print('wrote %s  (L=%.3f from mesh, %d interface faces)' % (out, L, len(tris)))


if __name__ == '__main__':
    main()
