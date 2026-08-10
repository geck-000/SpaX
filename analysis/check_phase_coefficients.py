# -*- coding: utf-8 -*-
"""Is the fitted brine coefficient really larger than the gas one?

The two are secants over very different ranges (brine 0.016-0.045, gas
0-0.099), and E/E_m is convex, so the comparison as fitted is not like for
like. Refit the gas sweep restricted to the brine range and see whether the
ordering survives.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')


def prep(f, keep=None):
    d = pd.read_csv(f)
    for c in ('E_x', 'E_matrix', 'phi_inclusion', 'porosity'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d = d.dropna(subset=['E_x', 'E_matrix'])
    if keep is not None:
        d = d[keep(d)]
    d['Erel'] = d.E_x / d.E_matrix
    return d


col = prep('results_colseeds_all.csv',
           keep=lambda d: d.run_id.str.extract(r'z(\d\d)')[0].astype(int) <= 75)
gas_all = prep('results_gas.csv')


def fit(d, label):
    X = np.column_stack([np.ones(len(d)), d.phi_inclusion, d.porosity])
    y = d.Erel.values
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    s2 = (r ** 2).sum() / (len(d) - 3)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    print('%-34s n=%3d  brine %.2f+-%.2f   gas %.2f+-%.2f'
          % (label, len(d), -b[1], se[1], -b[2], se[2]))
    return -b[1], -b[2], se[1], se[2]


print('as published (full gas range):')
both = pd.concat([col, gas_all], ignore_index=True)
b1, g1, sb1, sg1 = fit(both, '  pooled, gas to phi=0.099')

print()
print('gas restricted to the brine range (phi_gas <= 0.045):')
gas_m = gas_all[gas_all.porosity <= 0.045]
both_m = pd.concat([col, gas_m], ignore_index=True)
b2, g2, sb2, sg2 = fit(both_m, '  pooled, gas to phi=0.045')

print()
d1, sd1 = g1 - b1, np.hypot(sb1, sg1)
d2, sd2 = g2 - b2, np.hypot(sb2, sg2)
print('gas minus brine coefficient:')
print('   full gas range   : %+.3f +- %.3f  (%.1f sigma)' % (d1, sd1, abs(d1) / sd1))
print('   matched range    : %+.3f +- %.3f  (%.1f sigma)' % (d2, sd2, abs(d2) / sd2))
print()
print('Mori-Tanaka tangents put gas ABOVE brine by 0.16 (2.00 vs 1.84).')
print('A void cannot be less compliant than a filled pocket, so any fit that')
print('orders them the other way is reporting the fitting range, not the physics.')
