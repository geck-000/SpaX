#!/usr/bin/env python3
"""How much of the RVE-versus-field modulus offset is a gradient artefact?

Section 4.3.2 reports the homogenised moduli as ~2.3x stiffer than the
vibrating-beam fit of Marchenko (2024), and attributes the gap to a measurement
difference -- frequency dispersion and meso-scale compliance absent from the
periodic cell -- absorbed empirically by a 0.49 matrix factor.

Kujala et al. (1990) offer a second explanation for a discrepancy of that size.
A flexural test on a graded sheet does not return the local modulus; it returns
a rigidity-weighted average over a profile that is stiff at the top and soft at
the base. From a linear E(z) with alpha = E_bot/E_top they derive

    E_top / E_flex = 3(1+alpha) / (alpha^2 + 4 alpha + 1)                (K7)

which rises from 1 at alpha=1 to 2.3 at the alpha their beams exhibit. That is
numerically our offset, which is why the two explanations have to be separated
rather than assumed.

This script does the separation. It does not use the linear closed form, since
our own E(z) is convex (Section 4.3.3) and the closed form would misstate it.
Instead the same rigidity integrals are evaluated directly on the computed
profile:

    zbar = SUM E_i t_i z_i / SUM E_i t_i
    D    = SUM E_i [ t_i^3/12 + t_i (z_i - zbar)^2 ]
    E_flex = 12 D / H^3

so E_top/E_flex is whatever our profile actually implies. The question then is
how much of the 2.3x that accounts for, and how much is left for dispersion.

    python3 gradient_correction.py [column.csv ...]      # run from results/
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

OFFSET = 2.3          # reported RVE / vibrating-beam ratio, Section 4.3.2
BEAM_FACTOR = 0.49    # matrix factor currently used to absorb it


def load(path):
    """Mean E_eff per depth slice, ordered top to base."""
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        try:
            v = float(r['E_eff'])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[r['run_id'].split('_s')[0]].append(v)
    key = lambda k: int(''.join(c for c in k if c.isdigit()) or 0)
    return [st.mean(g[k]) for k in sorted(g, key=key)]


def flexural(E, H=1.0):
    """Rigidity-weighted flexural modulus of a stack of equal laminae."""
    E = np.asarray(E, float)
    n = len(E)
    t = H / n
    z = np.array([(i + 0.5) * t for i in range(n)])
    zbar = float((E * t * z).sum() / (E * t).sum())
    D = float((E * (t ** 3 / 12.0 + t * (z - zbar) ** 2)).sum())
    return 12.0 * D / H ** 3, zbar / H


def k7(alpha):
    """Kujala et al. (1990) Eq. (7): the linear-profile closed form."""
    return 3.0 * (1.0 + alpha) / (alpha ** 2 + 4.0 * alpha + 1.0)


def main():
    # The re-centred column is the production one: its base is the five-packing
    # mean, not the single reference packing that came out ~6 sigma low. Section
    # 4.3.3 assembles that column, so the gradient correction has to use it too
    # or the two sections disagree on alpha.
    paths = sys.argv[1:] or ['results_column_recentred.csv',
                             'results_steep_column.csv']
    print('%-28s %8s %9s %9s %9s %9s'
          % ('profile', 'alpha', 'E_top', 'E_flex', 'ratio', 'K7(alpha)'))
    rows = []
    for p in paths:
        if not os.path.isfile(p):
            continue
        E = load(p)
        if len(E) < 4:
            continue
        Ef, zb = flexural(E)
        a = E[-1] / E[0]
        rows.append((p, a, E[0], Ef, E[0] / Ef, k7(a), zb))
        print('%-28s %8.3f %9.3f %9.3f %9.3f %9.3f'
              % (os.path.basename(p), a, E[0] / 1e9, Ef / 1e9, E[0] / Ef, k7(a)))

    if not rows:
        raise SystemExit('no usable column files')

    print()
    print('The closed form K7 assumes a LINEAR E(z). Where it exceeds the ratio')
    print('computed from the actual profile, the profile is convex: it stays')
    print('stiff through the interior, so it is less gradient-affected than a')
    print('straight line between the same endpoints would be.')
    print()
    for p, a, Et, Ef, ratio, k, zb in rows:
        residual = OFFSET / ratio
        print('%s' % os.path.basename(p))
        print('   gradient alone explains        x%.2f  of the x%.1f offset'
              % (ratio, OFFSET))
        print('   residual for dispersion etc.   x%.2f' % residual)
        print('   that residual as a matrix factor: %.2f  (vs the %.2f in use)'
              % (1.0 / residual, BEAM_FACTOR))
        print('   neutral plane z0/H = %.3f' % zb)
        print()
    print('Interpretation. A residual close to x1 would mean the offset is')
    print('essentially all measurement geometry and the matrix factor is')
    print('double-counting. A residual close to x2.3 would mean the gradient')
    print('contributes nothing and the factor is doing real work.')


if __name__ == '__main__':
    main()
