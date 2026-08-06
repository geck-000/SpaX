"""Does the mesh contain the inclusion volume the deck asked for?

This is the gate on the paper re-run. The mesher used to build each inclusion
as its bounding sphere instead of the ellipsoid the packer placed, which
inflated the meshed inclusion fraction by 1/sphericity^2 -- 1.6x in the cold
cells, 2.2x in the warm ones. The moduli were correct for the cells built, but
the cells were not the ones the temperature and salinity specify.

The check is deliberately independent of both code paths it audits: the meshed
volume is summed straight from the deck's tetrahedra, and the target is read
from the deck's own parameter row.

    python3 audit_volume.py params/rve_x.csv out_x [tol]

Exits non-zero if the mean meshed/target ratio exceeds `tol` (default 1.15), so
a chained re-run stops rather than spending a thousand solves on bad geometry.
"""
import csv
import glob
import os
import re
import sys

import numpy as np


def parse_deck(path):
    """Node coordinates, tet connectivity and element sets from an .inp."""
    nodes, elems, esets = {}, {}, {}
    mode, cur = None, None
    for line in open(path, encoding='utf8', errors='replace'):
        s = line.strip()
        if not s:
            continue
        if s.startswith('*'):
            low = s.lower()
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
            else:
                mode = None
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


def meshed_fraction(path, L):
    nodes, elems, esets = parse_deck(path)
    if not elems:
        return None
    n = np.array([nodes[k] for k in sorted(nodes)])
    idx = {k: i for i, k in enumerate(sorted(nodes))}
    vol = {}
    for e, nd in elems.items():
        a, b, c, d = (n[idx[x]] for x in nd)
        vol[e] = abs(np.dot(np.cross(b - a, c - a), d - a)) / 6.0
    for name in ('Sphere_Only', 'SPHERE_ONLY'):
        if name in esets:
            ids = [i for i in expand(esets[name]) if i in vol]
            return sum(vol[i] for i in ids) / (L ** 3)
    return None


def main():
    deck_csv, outdir = sys.argv[1], sys.argv[2]
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1.15

    rows = {r['run_id']: r for r in csv.DictReader(
        open(deck_csv, encoding='utf8', errors='replace'))}

    ratios = []
    print('%-24s %10s %10s %8s' % ('run_id', 'target', 'meshed', 'ratio'))
    for f in sorted(glob.glob(os.path.join(outdir, 'Job-*-utx.inp')))[:40]:
        rid = os.path.basename(f)[4:].rsplit('-', 1)[0]
        r = rows.get(rid)
        if r is None:
            continue
        L = float(r['L'])
        # Target brine = meshed soft phase only. Gas voids are not meshed.
        incl = float(r.get('VoF_incl_sphere', 0) or 0)
        chan = float(r.get('channel_vof_target', 0) or 0)
        target = incl + chan
        got = meshed_fraction(f, L)
        if got is None or target <= 0:
            continue
        ratios.append(got / target)
        print('%-24s %10.4f %10.4f %8.3f' % (rid, target, got, got / target))

    if not ratios:
        print('audit: no comparable decks found -- not gating on this campaign')
        return 0
    m = float(np.mean(ratios))
    print('\nmean meshed/target over %d decks : %.3f   (tolerance %.2f)'
          % (len(ratios), m, tol))
    if m > tol:
        print('FAIL: the mesh still carries more inclusion than the deck asks for.')
        return 1
    print('PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main())
