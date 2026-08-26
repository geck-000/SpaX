# -*- coding: utf-8 -*-
"""Joint fit of the normalised knockdown across the two campaigns that vary the
two phases independently.

Normalising by each cell's own matrix modulus is what makes the campaigns
commensurable: the column's E_ice(T) falls 2.2% from surface to base, and left
in, that matrix softening is attributed to the brine.
"""
import os
import numpy as np
import pandas as pd

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')


def prep(f, keep=None):
    d = pd.read_csv(f)
    for c in ('E_x', 'E_matrix', 'phi_inclusion', 'porosity'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d = d.dropna(subset=['E_x', 'E_matrix'])
    if keep is not None:
        d = d[keep(d)]
    d['Erel'] = d['E_x'] / d['E_matrix']
    return d[['run_id', 'Erel', 'phi_inclusion', 'porosity']]


col = prep('results_colseeds_all.csv',
           keep=lambda d: d.run_id.str.extract(r'z(\d\d)')[0].astype(int) <= 75)
gas = prep('results_gas.csv')
both = pd.concat([col, gas], ignore_index=True)

for lab, d in (('column only  ', col), ('gas only     ', gas), ('POOLED       ', both)):
    X = np.column_stack([np.ones(len(d)), d.phi_inclusion, d.porosity])
    y = d.Erel.values
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    r2 = 1 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum()
    # standard errors
    dof = len(d) - 3
    s2 = (resid ** 2).sum() / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    print('%s n=%3d  intercept %.4f+-%.4f   brine %.2f+-%.2f   gas %.2f+-%.2f   R2=%.4f'
          % (lab, len(d), b[0], se[0], -b[1], se[1], -b[2], se[2], r2))

print()
print('range covered:')
print('  brine  %.4f - %.4f' % (both.phi_inclusion.min(), both.phi_inclusion.max()))
print('  gas    %.4f - %.4f' % (both.porosity.min(), both.porosity.max()))
print('  predictor correlation %.3f' % np.corrcoef(both.phi_inclusion, both.porosity)[0, 1])

# a single-coefficient form: are the two phases distinguishable at all?
X1 = np.column_stack([np.ones(len(both)), both.phi_inclusion + both.porosity])
y = both.Erel.values
b1, *_ = np.linalg.lstsq(X1, y, rcond=None)
r21 = 1 - ((y - X1 @ b1) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print()
print('single soft-phase coefficient: E/E_m = %.4f %+.2f phi_soft   R2=%.4f'
      % (b1[0], b1[1], r21))
