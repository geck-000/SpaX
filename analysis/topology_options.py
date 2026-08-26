# -*- coding: utf-8 -*-
"""What does a topology fix actually have to deliver?

Established: any ellipsoidal inclusion in a connected matrix leaves E_x pinned
near the Hashin-Shtrikman upper bound, because the ice percolates around every
inclusion however thin. Measured E_x at the base is 1.27 GPa against 6.5 from
the cell. The question is what geometry severs the load path enough, and by how
much a controlling parameter has to move.
"""
import numpy as np

E_ice, nu = 9.37, 0.33
E_br = 9 * 2.2 * 0.44e-3 / (3 * 2.2 + 0.44e-3)   # ~0.0013 GPa
PHI = 0.179
E_TARGET = 1.27

print('target: E_x = %.2f GPa, i.e. E_x/E_ice = %.4f' % (E_TARGET, E_TARGET / E_ice))
print('cell currently delivers 6.5 GPa (ratio %.3f)' % (6.5 / E_ice))
print()

# ---- (A) cell-spanning brine layers, interrupted by ice bridges
# In x the cell is ice in series with layer planes; within a plane a fraction b
# is ice bridge and carries the load, the rest is brine and carries none.
print('(A) brine LAYERS spanning the cell, with ice bridges of area fraction b')
print('    1/E_x = (1-t)/E_ice + t/(b E_ice),  t = layer thickness fraction')
print()
print('%10s %12s %14s' % ('t', 'b needed', 'bridge spacing'))
for t in (0.05, 0.10, 0.179, 0.30):
    # solve (1-t) + t/b = E_ice/E_target
    b = t / (E_ice / E_TARGET - (1 - t))
    print('%10.3f %12.4f %14s' % (t, b, '%.0f%% of plane is brine' % (100 * (1 - b))))
print()
print('    so the layers must be nearly complete: a few percent of the plane')
print('    carrying ice. That IS the skeletal layer near the interface, and it')
print('    is a parameter with a physical meaning tomography can measure.')
print()

# ---- (B) how sensitive is it? the danger of a knife-edge parameter
print('(B) sensitivity of E_x to the bridge fraction, at t = phi = %.3f' % PHI)
print('%10s %12s' % ('b', 'E_x GPa'))
for b in (0.01, 0.02, 0.03, 0.05, 0.10, 0.25, 0.50, 1.0):
    Ex = E_ice / ((1 - PHI) + PHI / b)
    print('%10.3f %12.3f' % (b, Ex))
print()
print('    E_x falls through the measured value between b = 0.02 and 0.03, and')
print('    is back above 4 GPa by b = 0.10. The model would be extremely')
print('    sensitive to a parameter that micro-CT resolves poorly, which is a')
print('    real argument AGAINST tuning it and FOR reporting the bound.')
print()

# ---- (C) what the current cell can and cannot span
print('(C) bounds at phi = %.3f, for reference' % PHI)
voigt = (1 - PHI) * E_ice + PHI * E_br
reuss = 1.0 / ((1 - PHI) / E_ice + PHI / E_br)
print('    Voigt %.3f   cell 6.50   measured %.2f   Reuss %.4f GPa'
      % (voigt, E_TARGET, reuss))
print('    the measurement sits at %.1f%% of the way from Reuss to Voigt on a'
      % (100 * (np.log(E_TARGET) - np.log(reuss)) / (np.log(voigt) - np.log(reuss))))
print('    log scale, so it is genuinely intermediate -- not near either limit.')
