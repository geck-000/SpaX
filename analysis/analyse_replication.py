# -*- coding: utf-8 -*-
"""Analyse the replication lane.

Three campaigns previously carried one packing per condition and the manuscript
had to hedge their conclusions; a fourth samples the gas fraction finely in the
dilute range to settle the phase-coefficient inversion. This reduces all four
and reports each claim with the scatter it now has.

Run from results/:  python3 ../analysis/analyse_replication.py
"""
import os
import numpy as np
import pandas as pd
from scipy import stats


def num(d, *cols):
    for c in cols:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    return d


def load(f, need=('E_x',)):
    if not os.path.exists(f):
        print('  %s not present yet' % f)
        return None
    d = num(pd.read_csv(f), 'E_x', 'E_z', 'phi_inclusion', 'porosity',
            'E_matrix', 'SCF_p99', 'MCnorm_p99', 'channel_vof_target')
    return d.dropna(subset=[c for c in need if c in d.columns])


def by_condition(d, pat):
    """Strip the _sN replicate suffix to recover the condition label."""
    d = d.copy()
    d['cond'] = d.run_id.astype(str).str.replace(r'_s\d+$', '', regex=True)
    return d


# ------------------------------------------------------------------ 1. failure
print('=' * 72)
print('FAILURE / CRITICAL ZONE  (was 1 packing per slice)')
print('=' * 72)
d = load('results_failure_rep.csv', need=('SCF_p99',))
if d is not None:
    d = by_condition(d, None)
    g = d.groupby('cond')['SCF_p99']
    print('%-10s %3s %9s %9s' % ('slice', 'n', 'P99', 'sd'))
    for c, v in g:
        print('%-10s %3d %9.3f %9.3f' % (c, len(v), v.mean(), v.std(ddof=0)))
    m = g.mean()
    base, top = m.iloc[-1], m.iloc[0]
    sd_base = g.std(ddof=0).iloc[-1]
    print()
    print('  base %.3f +- %.3f against surface %.3f' % (base, sd_base, top))
    print('  base exceeds surface by %.1f sd of the base ensemble'
          % ((base - top) / sd_base if sd_base else float('nan')))

# ------------------------------------------------------------------ 2/3. aniso
for f, lab, old in (('results_perc_rep.csv', 'PERCOLATION (was 2.4 sd)', 2.4),
                    ('results_chan_rep.csv', 'CHANNEL GEOMETRY (was 3.4 sd)', 3.4)):
    print()
    print('=' * 72)
    print('%s' % lab)
    print('=' * 72)
    d = load(f, need=('E_x', 'E_z'))
    if d is None:
        continue
    d = by_condition(d, None)
    d['a'] = d.E_z / d.E_x
    g = d.groupby('cond')['a']
    print('%-16s %3s %9s %9s' % ('condition', 'n', 'Ez/Ex', 'sd'))
    for c, v in g:
        print('%-16s %3d %9.4f %9.4f' % (c, len(v), v.mean(), v.std(ddof=0)))
    means, sds = g.mean(), g.std(ddof=0)
    span = means.max() - means.min()
    pooled = float(np.sqrt((sds ** 2).mean()))
    print()
    print('  span across conditions %.4f, pooled within-condition sd %.4f -> %.1f sd'
          % (span, pooled, span / pooled if pooled else float('nan')))
    print('  (previously %.1f sd on single cells)' % old)

# ------------------------------------------------------------------ 4. gas
print()
print('=' * 72)
print('DILUTE GAS SWEEP  (settles the phase-coefficient inversion)')
print('=' * 72)
d = load('results_gas_dilute.csv', need=('E_x',))
if d is not None:
    d['Erel'] = d.E_x / d.E_matrix
    X = np.column_stack([np.ones(len(d)), d.porosity, d.phi_inclusion])
    b, *_ = np.linalg.lstsq(X, d.Erel.values, rcond=None)
    r = d.Erel.values - X @ b
    s2 = (r ** 2).sum() / (len(d) - 3)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    print('  phi_gas range %.4f - %.4f over %d cells'
          % (d.porosity.min(), d.porosity.max(), len(d)))
    print('  E/E_m = %.4f - %.2f phi_gas - %.2f phi_brine' % (b[0], -b[1], -b[2]))
    print('  gas coefficient %.2f +- %.2f' % (-b[1], se[1]))
    print()
    print('  Mori-Tanaka tangent for a void is 2.00; the column brine fit gives 1.68.')
    print('  A gas coefficient above the brine one restores the physical ordering.')
