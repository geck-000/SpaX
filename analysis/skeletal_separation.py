# -*- coding: utf-8 -*-
"""Does arrangement matter, at skeletal-layer soft fractions?

Two models on the same 27 cells, each controlling for the realised soft
fraction and asking whether the morphology label adds anything:

    E_x    ~ 1 + phi_soft + is_channel_dominated
    E_z/E_x ~ 1 + phi_soft + is_channel_dominated

If the separation asserted for the column survives out here, the morphology
term should be null in the first and significant in the second.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')

d = pd.read_csv('results_skeletal.csv')
for c in ('E_x', 'E_z', 'phi_soft_total'):
    d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.dropna(subset=['E_x', 'E_z'])
d['chan'] = (d.run_id.str.extract(r'SKEL_([cp])')[0] == 'c').astype(float)
d['aniso'] = d.E_z / d.E_x
d['Ex'] = d.E_x / 1e9


def ols(y, cols):
    X = np.column_stack([np.ones(len(d))] + [d[c].values for c in cols])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ b
    dof = len(d) - X.shape[1]
    s2 = (resid ** 2).sum() / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = b / se
    p = 2 * (1 - stats.t.cdf(np.abs(t), dof))
    r2 = 1 - (resid ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return b, se, p, r2


print('n = %d cells, phi_soft %.3f-%.3f' % (len(d), d.phi_soft_total.min(),
                                            d.phi_soft_total.max()))
print()
for label, y in (('STIFFNESS  E_x [GPa]', d.Ex.values),
                 ('ANISOTROPY E_z/E_x  ', d.aniso.values)):
    b, se, p, r2 = ols(y, ['phi_soft_total', 'chan'])
    print('%s   R2=%.3f' % (label, r2))
    print('    intercept        %+8.4f +- %.4f   p=%.2g' % (b[0], se[0], p[0]))
    print('    phi_soft         %+8.4f +- %.4f   p=%.2g' % (b[1], se[1], p[1]))
    print('    channel-dominated%+8.4f +- %.4f   p=%.2g' % (b[2], se[2], p[2]))
    print()

# effect size in context
b, se, p, r2 = ols(d.aniso.values, ['phi_soft_total', 'chan'])
print('At matched soft fraction the channel-dominated cell is %.3f more'
      % b[2])
print('anisotropic, i.e. %.0f%% of the way from the column base (1.13) to the'
      % (100 * b[2] / (1.13 - 1.0)))
print('isotropic point -- a large effect at this scale.')
print()
b2, se2, p2, r22 = ols(d.Ex.values, ['phi_soft_total', 'chan'])
print('The same term in the stiffness model is %+.3f +- %.3f GPa (p=%.2f),'
      % (b2[2], se2[2], p2[2]))
print('against a %.1f GPa span of E_x across the sweep.'
      % (d.Ex.max() - d.Ex.min()))
