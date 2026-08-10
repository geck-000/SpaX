# -*- coding: utf-8 -*-
"""Direct per-cell couple-stress length from the second-order solve.

The size-sweep route infers a length scale from how the apparent modulus varies
with cell size, and needs a matched phi=0 control to remove the cube-versus-plate
extraction bias. The second-order scheme gives a more direct route: each cell
returns its own bending rigidity D_rve (moment conjugate to the prescribed
curvature at RP_K) AND its own first-order moduli, so the classical part can be
subtracted within the same cell:

    D_classical = [E_eff/(1-nu^2)] * L^4/12
    l^2         = (D_rve - D_classical) / (G_eff * L^2)

No regression, no reference cell, no size sweep. A real couple-stress length
gives l^2 > 0 in every cell independently.
"""
import os
import numpy as np
import pandas as pd

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')

D = 0.08     # mean inclusion diameter


def load(f):
    d = pd.read_csv(f)
    for c in ('L', 'E_eff', 'nu_eff', 'D_rve'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    return d.dropna(subset=['L', 'E_eff', 'nu_eff', 'D_rve'])


def lengths(d, label):
    L, E, nu, Dr = d.L.values, d.E_eff.values, d.nu_eff.values, d.D_rve.values
    G = E / (2 * (1 + nu))                       # isotropic estimate
    Dc = (E / (1 - nu ** 2)) * L ** 4 / 12.0
    l2 = (Dr - Dc) / (G * L ** 2)
    print('=== %s  (n=%d) ===' % (label, len(d)))
    print('   D_rve/D_classical : mean %.4f   range %.4f - %.4f'
          % ((Dr / Dc).mean(), (Dr / Dc).min(), (Dr / Dc).max()))
    print('   l^2 [m^2]         : mean %+.3e   %d of %d cells positive'
          % (l2.mean(), (l2 > 0).sum(), len(l2)))
    pos = l2[l2 > 0]
    if len(pos):
        print('   where positive, l : mean %.4f m = %.2f d' % (np.sqrt(pos).mean(),
                                                               np.sqrt(pos).mean() / D))
    # a bound: 95% upper limit on l from the mean of l^2 across cells
    se = l2.std(ddof=1) / np.sqrt(len(l2))
    hi = l2.mean() + 1.96 * se
    print('   mean l^2 = %+.3e +- %.3e (se); 95%% upper limit %+.3e' % (l2.mean(), se, hi))
    if hi > 0:
        print('   ->  l < %.4f m = %.2f d' % (np.sqrt(hi), np.sqrt(hi) / D))
    else:
        print('   ->  the whole interval is negative: no couple-stress length at any size')
    print()
    return l2


er = load('results_eringen.csv')
er['Ld'] = (er.L / D).round(0)
l2 = lengths(er, 'channelled base RVE, all cells')

print('by cell size:')
er['l2'] = l2
for Ld, g in er.groupby('Ld'):
    print('   L/d=%2.0f  n=%d  mean l^2 %+.3e   %d positive'
          % (Ld, len(g), g.l2.mean(), (g.l2 > 0).sum()))

print()
eh = load('results_eringen_homog.csv')
if len(eh):
    print('and the phi=0 control, which can carry no material length at all:')
    lengths(eh, 'homogeneous control')
