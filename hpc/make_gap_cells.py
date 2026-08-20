#!/usr/bin/env python3
r"""Resolve the sharp feature between phi = 0.0933 and 0.0970, or expose it.

The on-curve RAMP cells fall 3.575 -> 2.482 GPa across that interval, a 31 per
cent drop for a 4.5 per cent change in b, which is a local exponent of 7.86
against 0.60 on the interval below it and 2.79 on the interval above. Both cells
are internally tight -- seed scatter 0.016 and 0.011 GPa -- so it is not noise,
and the exponential fitted to G(phi) above phi_c has its length scale
lam = 0.0042 essentially set by this one gap. A closure feature resting on one
interval is the same mistake phi_sat = 0.104 was, so it is measured before it is
published.

GAP -- three brine fractions inside the interval, 0.094, 0.095 and 0.096, built
    exactly as the RAMP cells were: b from Assur at the realised target, slab
    fraction backed out through the pocket offset, drained only, the drained
    element sizing. Two seeds each rather than one, and the second seed is
    diagnostic rather than statistical here: the bridges are placed at random
    within each plane, so if the feature is a discrete event in the arrangement
    of two bridges, seeds should disagree across it in a way they have not
    disagreed anywhere else in this campaign.

NBR -- the same brine fraction, phi = 0.095, with four bridges to a plane
    instead of two, two seeds. This is the direct test of the geometric
    explanation. Everything else is held: same phi, same b, so the same areal
    ice fraction of the plane, just divided into four patches instead of two. If
    the sharp feature is material, it should not care. If it is an artefact of
    two wide patches in a periodic cell -- which is also what killed every
    b = 0.45 cell in LAYERB, three attempts apart -- then four narrower ones
    should shift or remove it.

    python3 make_gap_cells.py <outdir>
"""
import csv
import os
import sys

from make_layer_decks import BASE, COLS, K, LM_MAX, LM_MIN, assur_b, thickness

POCKET_OFFSET = 0.019      # measured across the LCOL and RAMP decks
ELEM_ACROSS = 1.6          # drained convergence, from rve_layermesh
MIN_DIST = '0.005'         # every cell here is below the thickness that failed


def mesh_for(slab, b):
    return min(max(thickness(slab, b) / ELEM_ACROSS, LM_MIN), LM_MAX)


def row(rid, phi, b, n_bridges):
    slab = round(phi - POCKET_OFFSET, 4)
    r = dict(BASE)
    r['run_id'] = rid
    r['slab_vof'] = '%.4f' % slab
    r['bridge_fraction'] = '%.4f' % b
    r['n_bridges'] = str(n_bridges)
    r['K_inclusion'] = K['drn']
    r['L_mesh'] = '%.4f' % mesh_for(slab, b)
    r['min_distance'] = MIN_DIST
    return [r[c] for c in COLS], slab, mesh_for(slab, b)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, 'rve_gapcells.csv')

    info = []
    with open(p, 'w', newline='', encoding='utf8') as fh:
        w = csv.writer(fh)
        w.writerow(COLS)
        # Inside the gap, two bridges, as RAMP had.
        for phi in (0.094, 0.095, 0.096):
            b = assur_b(phi)
            for s in (1, 2):
                rr, slab, lm = row('GAP_p%03d_drn_s%d' % (round(phi * 1000), s),
                                   phi, b, 2)
                w.writerow(rr)
            info.append((phi, b, slab, lm, 2))
        # The geometric control: same phi and b, four bridges.
        phi = 0.095
        b = assur_b(phi)
        for s in (1, 2):
            rr, slab, lm = row('NBR_p095_n4_drn_s%d' % s, phi, b, 4)
            w.writerow(rr)
        info.append((phi, b, slab, lm, 4))

    print('wrote %s  (%d cells)' % (p, 2 * len(info)))
    print('  %-8s %-8s %-9s %-9s %-9s %s'
          % ('phi', 'b', 'slab_vof', 'L_mesh', 'bridges', 'elements'))
    for phi, b, slab, lm, nb in info:
        print('  %-8.3f %-8.4f %-9.4f %-9.5f %-9d %.0fk'
              % (phi, b, slab, lm, nb, (0.5 / lm) ** 3 / 1e3))
    print()
    print('The interval being resolved: RAMP gives E = 3.575 GPa at phi = 0.0933')
    print('and 2.482 at 0.0970. These three sit inside it.')
    print()
    print('What the deck can return:')
    print('  a smooth fall through the gap  -> the feature is real and lam is')
    print('                                    measured rather than inferred')
    print('  a jump between two of them     -> something discrete; compare the')
    print('                                    four-bridge cells at the same phi')
    print('  seeds disagreeing in the gap   -> the bridge ARRANGEMENT matters,')
    print('                                    which makes it geometry not material')
    return 0


if __name__ == '__main__':
    sys.exit(main())
