"""Are the inclusions in a mesh ellipsoids, or spheres?

The mesher used to build every inclusion as its bounding sphere, discarding the
semi-axes the packer had placed. Reading the code tells you which branch is
live; this measures what is actually in the deck, which is the thing that
matters.

Each inclusion is recovered as a connected component of the Sphere_Only
element set (tetrahedra sharing nodes), and its three principal extents are
compared. A sphere gives 1.00 on every axis. A prolate ellipsoid at sphericity
s gives roughly s : s : 1, so the mean short/long ratio should land near the
deck's sphericity_avg rather than near unity.

Channels are excluded: they are cylinders spanning the cell and would swamp the
statistics.

    python3 check_inclusion_shape.py <deck.inp> [expected_sphericity]
"""
import collections
import re
import sys

import numpy as np


def parse(path):
    nodes, elems, esets = {}, {}, {}
    mode, cur = None, None
    for line in open(path, encoding='utf8', errors='replace'):
        s = line.strip()
        if not s:
            continue
        if s.startswith('*'):
            low = s.lower()
            mode = None
            if low.startswith('*node') and 'output' not in low:
                mode = 'node'
            elif low.startswith('*element') and 'output' not in low:
                mode = 'elem'
            elif low.startswith('*elset'):
                mode = 'elset'
                m = re.search(r'elset\s*=\s*([A-Za-z0-9_\-]+)', s, re.I)
                cur = m.group(1) if m else None
                if cur:
                    esets.setdefault(cur, {'gen': 'generate' in low, 'vals': []})
            continue
        p = [x for x in s.split(',') if x.strip()]
        if mode == 'node' and len(p) >= 4:
            nodes[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
        elif mode == 'elem' and len(p) >= 5:
            elems[int(p[0])] = [int(x) for x in p[1:5]]
        elif mode == 'elset' and cur:
            esets[cur]['vals'].extend(int(x) for x in p)
    return nodes, elems, esets


def expand(es):
    if not es['gen']:
        return es['vals']
    out, v = [], es['vals']
    for i in range(0, len(v) - 2, 3):
        out.extend(range(v[i], v[i + 1] + 1, v[i + 2]))
    return out


class UF(object):
    def __init__(self):
        self.p = {}

    def find(self, a):
        self.p.setdefault(a, a)
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    path = sys.argv[1]
    expect = float(sys.argv[2]) if len(sys.argv) > 2 else None
    nodes, elems, esets = parse(path)

    name = next((n for n in ('Sphere_Only', 'SPHERE_ONLY') if n in esets), None)
    if name is None:
        print('no inclusion element set in %s' % path)
        return 1
    ids = [e for e in expand(esets[name]) if e in elems]
    print('deck            : %s' % path)
    print('inclusion elems : %d' % len(ids))

    # connected components: elements sharing a node belong to one inclusion
    node_to_el = collections.defaultdict(list)
    for e in ids:
        for n in elems[e]:
            node_to_el[n].append(e)
    uf = UF()
    for e in ids:
        uf.find(e)
    for n, els in node_to_el.items():
        for e in els[1:]:
            uf.union(els[0], e)
    comps = collections.defaultdict(list)
    for e in ids:
        comps[uf.find(e)].append(e)

    def tetvol(e):
        a, b, c, d = (np.array(nodes[n]) for n in elems[e])
        return abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0

    L = max(max(c) for c in nodes.values())
    ratios, kept, skipped, merged, clipped = [], 0, 0, 0, 0
    for _root, els in comps.items():
        ns = {n for e in els for n in elems[e]}
        P = np.array([nodes[n] for n in ns])
        lo, hi = P.min(axis=0), P.max(axis=0)
        ext = hi - lo
        # a channel spans the cell in z; an inclusion does not
        if ext.max() > 0.5 * L:
            skipped += 1
            continue
        if len(ns) < 10 or ext.min() <= 0:
            continue
        # Only single, whole inclusions can report a meaningful axis ratio.
        # Two others must be excluded or the measurement is meaningless:
        #   - inclusions that overlap merge into one component, and a merged
        #     cluster is elongated whatever its members' shapes are;
        #   - an inclusion cut by a periodic face leaves a partial body.
        # Both are caught by comparing the component's volume with the
        # ellipsoid inscribed in its own bounding box: a single whole ellipsoid
        # fills it, a cluster or a fragment does not.
        vol = sum(tetvol(e) for e in els)
        vfill = vol / (np.pi / 6.0 * ext[0] * ext[1] * ext[2])
        if min(lo) < 1e-9 or max(hi) > L - 1e-9:
            clipped += 1
            continue
        if vfill < 0.90:
            merged += 1
            continue
        ratios.append(sorted(ext)[0] / sorted(ext)[2])
        kept += 1

    if not ratios:
        print('no resolvable inclusions found')
        return 1
    r = np.array(ratios)
    print('inclusions      : %d single whole bodies measured' % kept)
    print('                  (%d channel-like, %d face-clipped, %d merged '
          'clusters -- all excluded)' % (skipped, clipped, merged))
    print('short/long axis : mean %.3f   median %.3f   range %.3f-%.3f'
          % (r.mean(), np.median(r), r.min(), r.max()))
    print()
    if expect is not None:
        print('deck sphericity : %.3f' % expect)
    sph = r.mean() > 0.93
    print('VERDICT         : %s' % (
        'SPHERES -- mesher is discarding the semi-axes' if sph
        else 'ELLIPSOIDS -- the packer\'s semi-axes reached the mesh'))
    return 1 if sph else 0


if __name__ == '__main__':
    sys.exit(main())
