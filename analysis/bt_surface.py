# -*- coding: utf-8 -*-
"""Read the (b, t) grid: is the bridge exponent constant once t is controlled?

Runs under `abaqus python` (the ODB API), so Python 2.7 dialect and no imports
the Abaqus interpreter does not ship.

The closure writes the bridge factor as g = b^n with n constant. LAYERB found n
running 0.93, 1.03, 0.85 across b = 0.10, 0.20, 0.30 -- but that sweep held the
slab brine fraction, and phi_slab = t(1-b) means holding it while raising b
thickens the layer by 29%. Two variables moved, so the result cannot be read as
a b-dependence.

This grid fixes t along each row and varies b across it. Within a row, t is
constant by construction, so a variation in n there is a variation in b alone.

    abaqus python analysis/bt_surface.py <workdir>
"""
import glob
import math
import os
import re
import sys

from odbAccess import openOdb

E_ICE = 9.37
DRAIN = 1.04
L = 0.50
DISP = 0.005


def e_from_odb(path):
    """Volume-averaged S11 and the axial strain, straight from the last frame.

    Deliberately not routed through extract_first_order: this script is meant to
    be readable on its own and to depend on as little as possible, since it is
    the thing that decides whether a published functional form survives.
    """
    odb = openOdb(path, readOnly=True)
    try:
        step = odb.steps[odb.steps.keys()[0]]
        frame = step.frames[-1]
        if abs(frame.frameValue - 1.0) > 1e-9:
            return None, None
        S = frame.fieldOutputs['S']
        V = frame.fieldOutputs['EVOL']
        vol = {}
        for v in V.values:
            vol[v.elementLabel] = v.data
        num = 0.0
        den = 0.0
        for v in S.values:
            w = vol.get(v.elementLabel, 0.0)
            num += v.data[0] * w
            den += w
        if den <= 0:
            return None, None
        return num / den, den
    finally:
        odb.close()


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else '.'
    eps = DISP / L
    rows = {}
    for p in sorted(glob.glob(os.path.join(work, 'Job-BT_*-utx.odb'))):
        m = re.search(r'BT_t(\d+)_b(\d+)_s(\d+)', os.path.basename(p))
        if not m:
            continue
        t = int(m.group(1)) / 1000.0
        b = int(m.group(2)) / 100.0
        sig, vol = e_from_odb(p)
        if sig is None:
            print('%-28s INCOMPLETE' % os.path.basename(p))
            continue
        E = sig / eps / 1e9
        # phi for the pocket law: the slab brine plus the pocket population the
        # deck carries. Taken from the deck's own construction rather than read
        # back, because this script deliberately avoids the postprocessor.
        phi = t * (1.0 - b) + 0.019
        ep = E_ICE * (1.0 - 1.65 * phi) / DRAIN
        n = math.log(E / ep) / math.log(b)
        rows.setdefault((t, b), []).append((E, n))

    if not rows:
        print('no BT ODBs in %s' % work)
        return 1

    print('%6s %6s %9s %9s %9s' % ('t', 'b', 'E GPa', 'n', 'seeds'))
    table = {}
    for (t, b) in sorted(rows):
        vals = rows[(t, b)]
        Em = sum(v[0] for v in vals) / len(vals)
        nm = sum(v[1] for v in vals) / len(vals)
        table[(t, b)] = nm
        print('%6.2f %6.2f %9.4f %9.4f %9d' % (t, b, Em, nm, len(vals)))

    print('')
    print('n across b, at FIXED t (this is the controlled comparison):')
    worst = 0.0
    for t in sorted(set(k[0] for k in table)):
        ns = [table[(t, b)] for b in sorted(set(k[1] for k in table))
              if (t, b) in table]
        if len(ns) < 2:
            continue
        spread = max(ns) - min(ns)
        worst = max(worst, spread)
        print('  t=%.2f : %s   spread %.4f'
              % (t, ' '.join('%.4f' % x for x in ns), spread))

    print('')
    print('largest spread in n at fixed t: %.4f' % worst)
    print('LAYERB, which moved t as well:   0.1792')
    if worst < 0.05:
        print('')
        print('VERDICT: n is constant along b once t is held. b^n survives, and')
        print('the LAYERB result was the path through the plane, not the form.')
    else:
        print('')
        print('VERDICT: n still varies along b at fixed t. The bridge factor')
        print('needs both arguments; g(phi) cannot be recovered from b alone.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
