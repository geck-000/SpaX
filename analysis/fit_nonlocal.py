"""Fit the bending size sweep to both nonclassical families (plain python3).

Srinivasa & Reddy (Appl. Mech. Rev. 69(3) 030802, 2017) classify nonclassical
continua so that the sign of the size effect is a property of the theory class,
not of the material:

  type III.A  gradient / rotation-gradient (couple stress, MCST)
              the gradient is an energy penalty on curvature -> STIFFENING
              E_app/E_inf = 1 + 12 l^2 / L^2
  type III.B  integral nonlocal, Eringen kernel
              -> SOFTENING; deflections rise, frequencies and buckling loads fall
              E_inf/E_app = 1 + (e0a)^2 / L^2

The two are formally dual (exchange the roles of stress and strain), so a single
size sweep does not choose between them: it fixes the sign, and the sign selects
the family. The manuscript fits III.A, obtains a negative slope, and concludes
l^2 < 0. This script fits both, and subtracts the homogeneous-cube baseline
first so that whatever the cube-versus-plate kinematics and the discretisation
contribute is not read as material behaviour.

    python3 fit_nonlocal.py <channelled.csv> [homogeneous.csv] [--d 0.08]

Defaults to results_si2nd.csv (the published three-size sweep) and
results_homog.csv. Run from results/.
"""
import argparse
import csv
import os
import statistics as st
from collections import defaultdict

import numpy as np

NU = 0.33


def read(path, key='E_bending'):
    """Mean and population s.d. of a column, grouped by cell size.

    Population s.d. (ddof=0) throughout: every scatter quoted in these papers is
    over the packings of an ensemble, not a sample of a larger population.
    """
    g = defaultdict(list)
    for r in csv.DictReader(open(path)):
        try:
            v = float(r[key])
        except (ValueError, KeyError, TypeError):
            continue
        if v <= 0.0:                      # failed solves are written as 0/MISSING
            continue
        g[float(r['L'])].append(v)
    return {L: (st.mean(v), st.pstdev(v), len(v)) for L, v in sorted(g.items())}


def linfit(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    s, i = np.polyfit(x, y, 1)
    yh = i + s * x
    ss = float(((y - yh) ** 2).sum())
    tt = float(((y - y.mean()) ** 2).sum())
    return s, i, (1 - ss / tt if tt > 0 else float('nan'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('channelled', nargs='?', default='results_si2nd.csv')
    ap.add_argument('homog', nargs='?', default='results_homog.csv')
    ap.add_argument('--d', type=float, default=0.08, help='mean inclusion diameter')
    a = ap.parse_args()

    ch = read(a.channelled)
    if not ch:
        raise SystemExit('no usable rows in %s' % a.channelled)

    hom = read(a.homog) if os.path.isfile(a.homog) else {}

    print('channelled: %s' % a.channelled)
    print('%8s %6s %4s %12s %10s' % ('L', 'L/d', 'n', 'E_bend GPa', 'sd'))
    for L, (m, s, n) in ch.items():
        print('%8.2f %6.1f %4d %12.4f %10.4f' % (L, L / a.d, n, m / 1e9, s / 1e9))

    # Reference: the large-cell limit. Using the widest cell rather than a
    # separately computed plate modulus keeps the fit reference-free, which is
    # what the manuscript's slope method does.
    Ls = sorted(ch)
    E_ref = ch[Ls[-1]][0]
    r = {L: ch[L][0] / E_ref for L in Ls}

    if hom:
        print('\nhomogeneous baseline: %s' % a.homog)
        hL = sorted(hom)
        h_ref = hom[hL[-1]][0]
        print('%8s %12s %10s' % ('L', 'E_bend GPa', 'rel'))
        for L in hL:
            print('%8.2f %12.4f %10.4f' % (L, hom[L][0] / 1e9, hom[L][0] / h_ref))
        matched = [L for L in Ls if L in hom]
        if matched:
            print('  matched sizes: %s -> baseline subtracted point-by-point'
                  % ', '.join('%.2f' % L for L in matched))
            for L in matched:
                r[L] = r[L] / (hom[L][0] / h_ref)
        else:
            span = (hom[hL[-1]][0] / hom[hL[0]][0] - 1) * 100
            print('  NO matched sizes; baseline spans %+.1f%% over L=%.2f-%.2f'
                  % (span, hL[0], hL[-1]))
            print('  -> reported uncorrected; run rve_eringen_homog.csv for '
                  'matched-size subtraction')

    x = np.array([1.0 / L ** 2 for L in Ls])
    y_a = np.array([r[L] for L in Ls])           # III.A:  r     = 1 + 12 l^2/L^2
    y_b = np.array([1.0 / r[L] for L in Ls])     # III.B:  1/r   = 1 + (e0a)^2/L^2

    sA, iA, rA = linfit(x, y_a)
    sB, iB, rB = linfit(x, y_b)

    print('\n--- type III.A  (gradient / couple stress: stiffening) ---')
    print('  r = %.4f %+.4f / L^2      R2=%.4f' % (iA, sA, rA))
    if sA > 0:
        l2 = sA / 12.0
        print('  l = %.4f = %.2f d   (positive size effect)' % (l2 ** 0.5, l2 ** 0.5 / a.d))
    else:
        print('  slope < 0  ->  l^2 = %.4g < 0, no couple-stress length scale.' % (sA / 12.0))
        print('  This family is RULED OUT; it predicts the opposite sign.')

    print('\n--- type III.B  (integral nonlocal, Eringen: softening) ---')
    print('  1/r = %.4f %+.4f / L^2    R2=%.4f' % (iB, sB, rB))
    if sB > 0:
        e0a = sB ** 0.5
        print('  e0a = %.4f = %.2f d   (nonlocal length ~ %.1f%% of an inclusion diameter)'
              % (e0a, e0a / a.d, 100 * e0a / a.d))
    else:
        print('  slope < 0  ->  no softening nonlocality either.')
    print('  intercept %.4f (exact model gives 1.0000; departure measures how far'
          % iB)
    print('  the widest cell still is from the asymptote)')

    print('\nnote: with %d sizes spanning L/d = %.1f-%.1f, the intercept is'
          % (len(Ls), Ls[0] / a.d, Ls[-1] / a.d))
    print('constrained almost entirely by the widest cells. Treat a fitted')
    print('length as indicative until the sweep reaches the flat part of 1/L^2.')


if __name__ == '__main__':
    main()
