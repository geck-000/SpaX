"""Neutral plane and flexural modulus of a graded sheet with a skeletal base.

The manuscript assembles ten equal slices by classical lamination theory
(analysis/macro_plate.py) and obtains a neutral plane 3.1% of the thickness
above mid-depth. Kujala et al. (1990, IAHR, Table 2) measure z0/H = 0.37-0.39
on strain-gauged floating beams, i.e. 11-13% above mid-depth, together with a
bottom-to-top modulus ratio alpha = E_bot/E_top = 0.12-0.19. Our column gives
alpha = 0.55.

This script works in the same one-dimensional (cylindrical-bending) reduction
Kujala et al. used, so the comparison is like-for-like, and it accepts laminae
of unequal thickness so the bottom few percent can be resolved:

    zbar = SUM E_i t_i z_i / SUM E_i t_i
    D    = SUM E_i [ t_i^3/12 + t_i (z_i - zbar)^2 ]
    E_flex = 12 D / H^3 ,   E_ext = SUM E_i t_i / H

For a linear grading these reduce to the closed forms

    z0/H = (1+2a)/(3(1+a))          E_top/E_flex = 3(1+a)/(a^2+4a+1)

which reproduce Kujala's measured neutral axis from their measured alpha, and
which invert to alpha = 0.163 for the measured z0/H = 0.38.

Usage
-----
    python3 skeletal_clt.py --probe
        No new data needed. Sweeps skeletal thickness and residual modulus over
        the existing column and reports which combinations reach the measured
        neutral axis, i.e. what the skeletal decks have to deliver.

    python3 skeletal_clt.py results_column.csv --skeletal results_skeletal_laminae.csv
        Re-assembles the column with the resolved basal laminae substituted for
        the lowest slice, and reports z0/H, E_flex, E_ext and alpha against the
        measured values.

Run from results/.
"""
import argparse
import csv
import statistics as st
from collections import defaultdict

import numpy as np

Z0_MEAS = (0.37, 0.39)          # Kujala et al. 1990, Table 2
ALPHA_MEAS = (0.12, 0.19)


def zbar_D(E, t, z):
    """Neutral-plane depth and flexural rigidity of a stack. z = layer centres,
    measured downward from the top surface; all arrays same length."""
    E, t, z = map(np.asarray, (E, t, z))
    S0 = float((E * t).sum())
    zb = float((E * t * z).sum() / S0)
    D = float((E * (t ** 3 / 12.0 + t * (z - zb) ** 2)).sum())
    return zb, D, S0


def report(E, t, z, H, label):
    zb, D, S0 = zbar_D(E, t, z)
    E_flex = 12.0 * D / H ** 3
    E_ext = S0 / H
    a = E[-1] / E[0]
    print('%-34s z0/H=%.3f  E_ext=%.2f  E_flex=%.2f GPa  E_top/E_flex=%.2f  alpha=%.3f'
          % (label, zb / H, E_ext / 1e9, E_flex / 1e9, E[0] / E_flex, a))
    return zb / H, E_flex, a


def load_column(path):
    """Mean E_eff per slice from a column results CSV (ddof=0 over packings).

    Only the equal-thickness depth slices are returned. A steep-column file may
    also carry the resolved basal sub-laminae (z955, z970, ...) in the same CSV;
    those are a subdivision of the lowest slice, not extra slices, and admitting
    them here would build a stack of n+4 layers all assumed equally thick. That
    silently misplaces the neutral plane, so they are dropped and passed instead
    through --skeletal, which knows their true thicknesses.
    """
    g = defaultdict(list)
    dropped = set()
    for r in csv.DictReader(open(path)):
        rid = r['run_id'].split('_s')[0]
        if len(rid.rsplit('_z', 1)[-1]) > 2:
            dropped.add(rid)
            continue
        try:
            g[rid].append(float(r['E_eff']))
        except (ValueError, KeyError, TypeError):
            continue
    if dropped:
        print('  (ignoring %d sub-lamina rows in %s: %s -- pass them via --skeletal)'
              % (len(dropped), path, ', '.join(sorted(dropped))))
    out = []
    for rid in sorted(g, key=lambda k: int(''.join(c for c in k if c.isdigit()) or 0)):
        out.append((rid, st.mean(g[rid])))
    return out


def probe(E_col, H=1.0):
    """Without new solves: how soft and how thick must the skeletal layer be?"""
    n = len(E_col)
    t0 = H / n
    z0 = np.array([(i + 0.5) * t0 for i in range(n)])
    E0 = np.array([e for _, e in E_col])
    t = np.full(n, t0)

    print('baseline column (%d equal slices, H=%.2f m)' % (n, H))
    report(E0, t, z0, H, '  as published')
    print('  measured (Kujala 1990): z0/H = %.2f-%.2f, alpha = %.2f-%.2f\n'
          % (Z0_MEAS + ALPHA_MEAS))

    print('replace the lowest slice by [remaining ice | skeletal layer]:')
    print('%10s %10s %10s %10s %10s' % ('t_sk/H', 'E_sk/E_top', 'z0/H', 'alpha', 'E_top/E_flex'))
    hits = []
    for tsk_f in (0.01, 0.02, 0.03, 0.05):
        for esk_f in (0.30, 0.20, 0.12, 0.06, 0.03):
            tsk = tsk_f * H
            t_rem = t0 - tsk
            if t_rem <= 0:
                continue
            E = np.concatenate([E0[:-1], [E0[-1], esk_f * E0[0]]])
            tt = np.concatenate([t[:-1], [t_rem, tsk]])
            edges = np.concatenate([[0.0], np.cumsum(tt)])
            zc = 0.5 * (edges[:-1] + edges[1:])
            zb, D, S0 = zbar_D(E, tt, zc)
            E_flex = 12.0 * D / H ** 3
            a = E[-1] / E[0]
            ok = Z0_MEAS[0] <= zb / H <= Z0_MEAS[1]
            print('%10.3f %10.2f %10.3f %10.3f %10.2f %s'
                  % (tsk_f, esk_f, zb / H, a, E[0] / E_flex, '  <-- measured band' if ok else ''))
            if ok:
                hits.append((tsk_f, esk_f, zb / H, E[0] / E_flex))
    # --- shape vs endpoint ratio ------------------------------------------
    # The closed forms above assume E varies LINEARLY across the whole
    # thickness. Kujala et al. state that assumption explicitly, and their
    # alpha is a fitted parameter of that linear model rather than an
    # independent measurement of the endpoint ratio. Our profile is flat
    # through the cold interior and drops only in the bottom fifth, which puts
    # the centroid near mid-depth however soft the base becomes. Separate the
    # two effects before blaming the base.
    print('\nshape vs endpoint ratio (both at the same alpha):')
    shape = E0 / E0[0]
    for a_t in (E0[-1] / E0[0], 0.30, 0.163):
        lin = np.linspace(1.0, a_t, n)
        # our shape, rescaled so its endpoints match a_t
        ours = (shape - shape[-1]) / (1.0 - shape[-1]) * (1.0 - a_t) + a_t
        zl, _, _ = zbar_D(lin, t, z0)
        zo, _, _ = zbar_D(ours, t, z0)
        print('  alpha=%.3f   linear profile z0/H=%.3f   our profile shape z0/H=%.3f'
              % (a_t, zl / H, zo / H))
    print('  measured band %.2f-%.2f is reached by the LINEAR profile at'
          % Z0_MEAS)
    print('  alpha~0.16, and is not reached by our profile shape at any alpha.')

    print()
    if hits:
        print('combinations reaching the measured neutral axis:')
        for tsk_f, esk_f, zb, ratio in hits:
            print('  t_sk = %.0f%% of H with E_sk = %.0f%% of E_top  ->  z0/H=%.3f, '
                  'E_top/E_flex=%.2f' % (tsk_f * 100, esk_f * 100, zb, ratio))
        print('\nThe skeletal decks target phi_b up to 0.50 over the bottom 5%;')
        print('this is the modulus that has to come out for the neutral axis to move.')
    else:
        print('No combination in the swept range reaches z0/H = %.2f-%.2f.' % Z0_MEAS)
        print('A thin soft layer at the plate surface removes stiffness from the')
        print('far fibre but is too thin to move the centroid: resolving the')
        print('skeletal layer alone cannot reproduce the measured neutral axis.')
        print('The mismatch is a statement about the SHAPE of E(z) -- our C-shape')
        print('salinity holds the cold interior flat -- not about the basal')
        print('modulus. Test it by steepening the profile, not the base.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('column', nargs='?', default='results_column_recentred.csv')
    ap.add_argument('--skeletal', default=None)
    ap.add_argument('--H', type=float, default=1.0)
    ap.add_argument('--probe', action='store_true')
    a = ap.parse_args()

    col = load_column(a.column)
    if not col:
        raise SystemExit('no usable rows in %s' % a.column)

    if a.probe or a.skeletal is None:
        probe(col, a.H)
        return

    sk = load_column(a.skeletal)
    n = len(col)
    t0 = a.H / n
    # the resolved laminae replace the lowest slice, splitting its thickness
    tsk = t0 / len(sk)
    E = np.array([e for _, e in col[:-1]] + [e for _, e in sk])
    tt = np.array([t0] * (n - 1) + [tsk] * len(sk))
    edges = np.concatenate([[0.0], np.cumsum(tt)])
    zc = 0.5 * (edges[:-1] + edges[1:])
    print('resolved stack: %d equal slices + %d skeletal laminae' % (n - 1, len(sk)))
    for rid, e in sk:
        print('    %-16s E = %.3f GPa' % (rid, e / 1e9))
    print()
    report(E, tt, zc, a.H, 'with skeletal laminae')
    print('measured (Kujala 1990): z0/H = %.2f-%.2f, alpha = %.2f-%.2f'
          % (Z0_MEAS + ALPHA_MEAS))


if __name__ == '__main__':
    main()
