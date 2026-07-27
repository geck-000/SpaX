#!/usr/bin/env python3
"""Compare the warm-base full-tensor ensembles across cell size.

Two five-packing ensembles of the same base slice, solved in all six load
cases, differing only in the cell edge: L=0.50 (3-5 channels per cell) and
L=0.80 (10-11). The comparison separates two things that a single cell size
cannot:

  * ratios -- E_y/E_x, E_z/E_xy, G_axial/G_xy. These are within-cell
    comparisons, so systematic errors largely cancel and they converge fast.
  * the modulus itself -- E_x. This is an absolute quantity and, at the base
    soft-phase fraction (~30% of the cell, close to the continuum percolation
    threshold for overlapping spheres), it is strongly size-sensitive.

The in-plane ratio is the one the ensembles were run to settle: at L=0.50 it
sat at 1.012 with all five packings on the same side of unity, which looked
systematic; at L=0.80 it is 0.998 with the packings split, i.e. in-plane
isotropy holds and the L=0.50 offset was a small-cell artifact.

Usage: python3 compare_basetensor_sizes.py [results_L050.csv results_L080.csv]
"""
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULTS = (os.path.join(HERE, '..', 'results', 'results_basetensor_seeds.csv'),
            os.path.join(HERE, '..', 'results', 'results_basetensor_bt80.csv'))
LABELS = ('L=0.50', 'L=0.80')


def stat(df, key, scale=1.0):
    v = df[key].values / scale
    return v.mean(), v.std(ddof=0)


def welch(a, b):
    """Separation of two means in sigma, Welch style (both n=5)."""
    sa = a.std(ddof=1) ** 2 / len(a)
    sb = b.std(ddof=1) ** 2 / len(b)
    d = abs(a.mean() - b.mean())
    return d / math.sqrt(sa + sb) if (sa + sb) > 0 else float('inf')


def main():
    paths = sys.argv[1:3] if len(sys.argv) > 2 else DEFAULTS
    d = [pd.read_csv(p) for p in paths]
    if any(len(x) == 0 for x in d):
        print('empty input')
        return 1

    print('warm-base full-tensor ensembles, %d and %d packings\n'
          % (len(d[0]), len(d[1])))
    print('  %-24s %-22s %-22s %s'
          % ('quantity', LABELS[0], LABELS[1], 'separation'))
    rows = (('E_x', 'E_x (GPa)', 1e9), ('E_y', 'E_y (GPa)', 1e9),
            ('E_z', 'E_z (GPa)', 1e9),
            ('inplane_ratio', 'E_y/E_x', 1.0),
            ('E_ratio', 'E_z/E_xy', 1.0),
            ('G_ratio', 'G_axial/G_xy', 1.0))
    for key, lab, sc in rows:
        m0, s0 = stat(d[0], key, sc)
        m1, s1 = stat(d[1], key, sc)
        sep = welch(d[0][key].values / sc, d[1][key].values / sc)
        print('  %-24s %8.4f +/- %-9.4f %8.4f +/- %-9.4f %.1f sigma'
              % (lab, m0, s0, m1, s1, sep))

    print('\nin-plane isotropy (the question these were run to settle):')
    for lab, x in zip(LABELS, d):
        r = x.inplane_ratio.values
        sem = r.std(ddof=1) / math.sqrt(len(r))
        print('  %-7s E_y/E_x = %.4f, %d/%d packings above unity, '
              'mean %.1f SEM from it'
              % (lab, r.mean(), int((r > 1).sum()), len(r),
                 abs(r.mean() - 1) / sem if sem else float('inf')))

    e0, e1 = d[0].E_x.values.mean(), d[1].E_x.values.mean()
    print('\nabsolute modulus: E_x %.3f -> %.3f GPa (%+.1f%%)'
          % (e0 / 1e9, e1 / 1e9, 100 * (e1 - e0) / e0))
    print('  the ratios above are cell-size robust; E_x is not. At the base '
          'soft-phase\n  fraction the cell sits near percolation, where the '
          'modulus keeps falling with\n  cell size. The published box-size '
          'sweep (results_sizechan.csv) shows E_x flat\n  to 0.3% over the '
          'same L range, but at roughly a third of this soft fraction,\n  so '
          'it does not certify convergence here.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
