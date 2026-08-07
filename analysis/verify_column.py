#!/usr/bin/env python3
"""Recompute the depth-column chain and check it against the manuscript.

Companion to verify_sizeeffect.py, which covers the bending study. This one
covers everything that descends from the production column: Table 2, the
replicate scatter and anisotropy significance of Section 4.1.1, the softening
narrative of Section 4.1.2, the laminated plate of Section 4.3.3 and the
gradient separation of Section 4.3.2.

It exists because that chain is long and every link is quoted in the text. The
column feeds the CLT, which feeds the neutral plane and the bend-stretch
coupling; it also feeds the flexural ratio, which feeds the residual and hence
the matrix factor. Changing the column -- as switching from the single
reference packing to the five-packing ensemble mean did -- moves all of them at
once, and a number missed in propagation is invisible by inspection.

Every quantity is compared against the value in the manuscript and the script
exits non-zero if any has drifted.

    python3 verify_column.py            # run from results/
    python3 verify_column.py --verbose  # also print the per-slice table
"""
import argparse
import csv
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

COLUMN = 'results_column_ensemble.csv'
SEEDS = 'results_colseeds_all.csv'
BEAM_FACTOR = 0.49
# Column-mean ratio of the uncalibrated RVE column to the field profile.
# 2.3 before the Marchenko curve was corrected to his Eq. (17); 3.11 after.
OFFSET = 3.11

# Marchenko (2024) Eq. (17), p.11, as used by plot_marchenko_match.py.
# E0 is the surface modulus and alpha the base/surface ratio.
M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6

# (label, quoted value, tolerance). Section numbers refer to the manuscript.
QUOTED = {
    # Table 2 and Section 4.1.2
    'E_surface':        (8.92, 0.005),
    'E_peak':           (9.06, 0.005),
    'E_base':           (6.50, 0.006),
    'peak_depth':       (0.25, 0.001),
    'surface_below_pct': (1.5, 0.05),
    'knockdown_pct':    (28.3, 0.05),
    # Section 4.1.1
    'cov_surface_pct':  (0.21, 0.005),
    'cov_base_pct':     (1.72, 0.005),
    'cov_max_pct':      (1.72,  0.05),
    'aniso_z85':        (1.0505, 0.0005),
    'aniso_z85_sigma':  (4.9,  0.05),
    'aniso_base':       (1.1324, 0.0005),
    'aniso_base_sigma': (10.2,  0.05),
    # Section 4.3.3, laminated plate
    'E_ext':            (8.58, 0.005),
    'E_flex':           (8.18, 0.005),
    'E_ext_beam':       (4.20, 0.005),
    'E_flex_beam':      (4.01, 0.005),
    'neutral_offset_pct': (1.9, 0.05),
    'coupling':         (0.068, 0.005),
    # Section 4.3.2, gradient separation
    'alpha':            (0.73, 0.005),
    'grad_ratio':       (1.10, 0.005),
    'k7_alpha':         (1.17, 0.005),
}


def load_column(path):
    rows = list(csv.DictReader(open(path)))
    rows.sort(key=lambda r: int(''.join(c for c in r['run_id'] if c.isdigit())))
    return rows


def q_plane_stress(Ex, Ey, nu_xy, Gxy):
    nu_yx = nu_xy * Ey / Ex
    d = 1.0 - nu_xy * nu_yx
    Q = np.zeros((3, 3))
    Q[0, 0] = Ex / d
    Q[1, 1] = Ey / d
    Q[0, 1] = Q[1, 0] = nu_xy * Ey / d
    Q[2, 2] = Gxy
    return Q


def clt(rows, H=1.0, scale=1.0):
    """A, B, D of the stack; returns the engineering quantities the paper quotes."""
    n = len(rows)
    zk = np.linspace(-H / 2.0, H / 2.0, n + 1)
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    for k, r in enumerate(rows):
        Q = q_plane_stress(float(r['E_x']) * scale, float(r['E_y']) * scale,
                           float(r.get('nu_x') or 0.33), float(r['G_xy']) * scale)
        z1, z2 = zk[k], zk[k + 1]
        A += Q * (z2 - z1)
        B += Q * 0.5 * (z2 ** 2 - z1 ** 2)
        D += Q * (1.0 / 3.0) * (z2 ** 3 - z1 ** 3)
    E_ext = 1.0 / (np.linalg.inv(A)[0, 0] * H)
    E_flex = 12.0 / (np.linalg.inv(D)[0, 0] * H ** 3)
    z_na = B[0, 0] / A[0, 0]
    coupling = abs(B[0, 0]) / np.sqrt(A[0, 0] * D[0, 0])
    return E_ext, E_flex, z_na / H, coupling


def rigidity(E, H=1.0):
    """Neutral plane and flexural modulus of equal laminae (the 1-D reduction
    Kujala et al. use, so the comparison is like for like)."""
    E = np.asarray(E, float); n = len(E); t = H / n
    z = np.array([(i + 0.5) * t for i in range(n)])
    zb = float((E * t * z).sum() / (E * t).sum())
    D = float((E * (t ** 3 / 12.0 + t * (z - zb) ** 2)).sum())
    return zb / H, 12.0 * D / H ** 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()
    for f in (COLUMN, SEEDS):
        if not os.path.isfile(f):
            raise SystemExit('missing %s -- run from results/' % f)

    rows = load_column(COLUMN)
    E = np.array([float(r['E_x']) for r in rows]) / 1e9
    z = np.arange(0.05, 1.0, 0.1)
    got = {}

    # ---- column shape ----------------------------------------------------
    pk = int(np.argmax(E))
    got['E_surface'] = E[0]
    got['E_peak'] = E[pk]
    got['E_base'] = E[-1]
    got['peak_depth'] = z[pk]
    got['surface_below_pct'] = 100.0 * (1.0 - E[0] / E[pk])
    got['knockdown_pct'] = 100.0 * (1.0 - E[-1] / E[pk])

    drops = np.diff(E[pk:])
    mono = bool((drops < 0).all())
    steep = bool((np.diff(np.abs(drops)) > 0).all())

    # ---- replicate scatter and anisotropy --------------------------------
    g = defaultdict(list)
    for r in csv.DictReader(open(SEEDS)):
        g[r['run_id'].split('_s')[0]].append(r)
    keys = sorted(g, key=lambda k: int(''.join(c for c in k if c.isdigit())))
    covs, ratios = [], {}
    for k in keys:
        ex = [float(x['E_x']) for x in g[k]]
        rr = [float(x['E_z']) / float(x['E_x']) for x in g[k]]
        covs.append(100.0 * st.pstdev(ex) / st.mean(ex))
        ratios[k] = (st.mean(rr), st.pstdev(rr))
    got['cov_surface_pct'] = covs[0]
    got['cov_base_pct'] = covs[-1]
    got['cov_max_pct'] = max(covs)
    n_under2 = sum(1 for c in covs if c < 2.0)
    for tag, key in (('z85', keys[-2]), ('base', keys[-1])):
        m, s = ratios[key]
        got['aniso_%s' % tag] = m
        got['aniso_%s_sigma' % tag] = (m - 1.0) / s

    # ---- laminated plate --------------------------------------------------
    E_ext, E_flex, z_na, coup = clt(rows)
    got['E_ext'], got['E_flex'] = E_ext / 1e9, E_flex / 1e9
    got['neutral_offset_pct'] = abs(z_na) * 100.0
    got['coupling'] = coup
    Eb_ext, Eb_flex, _, _ = clt(rows, scale=BEAM_FACTOR)
    got['E_ext_beam'], got['E_flex_beam'] = Eb_ext / 1e9, Eb_flex / 1e9

    # ---- gradient separation ---------------------------------------------
    Eeff = np.array([float(r['E_eff']) for r in rows]) / 1e9
    zb, Ef = rigidity(Eeff)
    a = Eeff[-1] / Eeff[0]
    got['alpha'] = a
    got['grad_ratio'] = Eeff[0] / Ef
    got['k7_alpha'] = 3.0 * (1 + a) / (a ** 2 + 4 * a + 1)
    got['residual'] = OFFSET / (Eeff[0] / Ef)
    got['matrix_factor'] = 1.0 / got['residual']


    if args.verbose:
        print('%6s %9s %9s %8s %9s' % ('z/H', 'E_x', 'CoV%', 'Ez/Ex', 'sigma'))
        for i, k in enumerate(keys):
            m, s = ratios[k]
            print('%6.2f %9.3f %9.3f %8.4f %9.2f'
                  % (z[i], E[i], covs[i], m, (m - 1) / s if s else 0))
        print()

    print('%-22s %11s %11s   %s' % ('quantity', 'computed', 'in the text', ''))
    bad = []
    for k, (want, tol) in QUOTED.items():
        have = got[k]
        ok = abs(have - want) <= tol
        if not ok:
            bad.append((k, have, want))
        print('%-22s %11.4f %11.4f   %s'
              % (k, have, want, 'ok' if ok else '<-- DIFFERS'))

    print()
    print('Narrative claims that are not single numbers:')
    print('  below the peak the softening is monotonic          : %s' % mono)
    print('  ... and monotonically steepening                   : %s' % steep)
    print('  all ten slices below 2%% CoV                       : %s (%d)'
          % (n_under2 == 10, n_under2))
    if not (mono and steep and n_under2 == 10):
        bad.append(('narrative claim', 0, 0))

    if bad:
        print('\n%d quantity(ies) no longer match the manuscript:' % len(bad))
        for k, have, want in bad:
            print('   %-22s computed %.4f, text says %.4f' % (k, have, want))
        return 1
    print('\nAll %d quoted quantities reproduce, and the three narrative claims hold.'
          % len(QUOTED))
    return 0


if __name__ == '__main__':
    sys.exit(main())
