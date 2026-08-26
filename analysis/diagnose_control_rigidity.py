# -*- coding: utf-8 -*-
"""Why does the phi=0 control return D_rve > D_classical?

For the control every input to Eq. (direct) is exact: the cell is genuinely
homogeneous and isotropic, so mu = E/2(1+nu) is not an approximation and
E_eff, nu_eff are the matrix values. Any excess therefore comes from the
boundary value problem, not from the arithmetic. Candidates:

  (a) the wrong classical reference (plate vs 3D periodic constraint)
  (b) the Lesicar integral constraints, which are kinematic and can only stiffen
  (c) discretisation
"""
import os
import numpy as np
import pandas as pd

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')

d = pd.read_csv('results_eringen_homog.csv')
for c in ('L', 'E_eff', 'nu_eff', 'D_rve', 'E_matrix'):
    d[c] = pd.to_numeric(d[c], errors='coerce')
d = d.dropna(subset=['L', 'E_eff', 'nu_eff', 'D_rve']).sort_values('L')

print('control cells: E_eff vs the matrix modulus they were built from')
print('%8s %12s %12s %10s' % ('L', 'E_eff GPa', 'E_matrix GPa', 'ratio'))
for _, r in d.iterrows():
    print('%8.2f %12.4f %12.4f %10.5f'
          % (r.L, r.E_eff / 1e9, r.E_matrix / 1e9, r.E_eff / r.E_matrix))
print()
print('-> first-order extraction on the control is exact to %.2g, so E_eff and'
      % abs(d.E_eff / d.E_matrix - 1).max())
print('   nu_eff are not the problem.')
print()

L, E, nu, Dr = d.L.values, d.E_eff.values, d.nu_eff.values, d.D_rve.values
I = L ** 4 / 12.0
print('classical reference, three candidate constraints:')
for lab, Dc in (('plane stress  E/(1-nu^2)', (E / (1 - nu ** 2)) * I),
                ('uniaxial      E',          E * I),
                ('plane strain  E(1-nu)/((1+nu)(1-2nu))',
                 E * (1 - nu) / ((1 + nu) * (1 - 2 * nu)) * I)):
    print('   %-38s D_rve/D_c = %.4f' % (lab, (Dr / Dc).mean()))
print()
print('-> the measured rigidity sits between the plane-stress and plane-strain')
print('   references, i.e. the periodic cube is stiffer in bending than a plate')
print('   of the same material but softer than one held in plane strain.')
print()

# does the excess depend on cell size? a kinematic-constraint artefact should
# be roughly size-independent; a discretisation one should shrink with L/L_mesh
d['ratio'] = Dr / ((E / (1 - nu ** 2)) * I)
print('excess against cell size (mesh size fixed at L_mesh=0.033):')
for _, r in d.iterrows():
    print('   L=%.2f  L/L_mesh=%4.1f   D_rve/D_classical = %.4f'
          % (r.L, r.L / 0.033, r.ratio))
sl = np.polyfit(d.L / 0.033, d.ratio, 1)[0]
print('   trend with elements-per-edge: %+.2e per element' % sl)
