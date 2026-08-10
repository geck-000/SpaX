# -*- coding: utf-8 -*-
"""Why the brine layers are too stiff: the fill is modelled as gas-free.

The bulk-modulus sweep settled the mechanism. At fixed geometry, dropping the
fill's bulk modulus from 2.2 GPa to 0.022 GPa takes E_x from 4.605 to 0.780 GPa
while E_z stays at 7.49 within 0.6%. A spanning layer is fully confined, so it
cannot thin under normal stretch without changing volume and therefore resists
at K, not at E. The ice bridges never controlled anything.

So the question is not what b should be but what K the fill really has. Brine
is modelled here as pure salt water, K = 2.2 GPa. Real brine is not gas-free:
dissolved air comes out of solution as ice freezes and collects with the brine
it is expelled alongside. A liquid carrying even a trace of free gas is far
more compressible than the liquid alone, because the gas takes nearly all of
the volume change. Wood's law gives the mixture,

    1/K_mix = f/K_gas + (1-f)/K_brine,

and with K_gas about 1.4e-4 GPa for air at atmospheric pressure the collapse is
violent: tenths of a percent of gas cost an order of magnitude in K.

This is not a free parameter dressed up. Sea ice gas content is measured, and
the cells already carry a gas fraction -- but as SEPARATE voids elsewhere in
the cell, which leaves the brine incompressible. Whether the gas sits apart
from the brine or within it is a topological statement about the microstructure
and it changes the modulus by a factor of several.
"""
import numpy as np

K_BRINE = 2.2            # GPa, salt water
K_GAS = 1.4e-4           # GPa, air at 1 atm, adiabatic (gamma*P)

# bulk-modulus sweep, this session: identical geometry, b = 0.03, phi = 0.179
SWEEP_K = np.array([2.2, 0.22, 0.022])
SWEEP_EX = np.array([4.605, 1.612, 0.780])
E_MEASURED = 1.27        # GPa, Kujala base, mean of four beams


def wood(f, k_gas=K_GAS, k_liq=K_BRINE):
    """Bulk modulus of a liquid carrying free-gas volume fraction f."""
    f = np.asarray(f, dtype=float)
    return 1.0 / (f / k_gas + (1.0 - f) / k_liq)


def K_for(E_target):
    """Interpolate the sweep, in log-log, for the K that gives E_target."""
    o = np.argsort(SWEEP_EX)
    return float(np.exp(np.interp(np.log(E_target),
                                  np.log(SWEEP_EX[o]), np.log(SWEEP_K[o]))))


def main():
    print('bulk-modulus sweep (b = 0.03, identical geometry)')
    print('%12s %12s' % ('K fill GPa', 'E_x GPa'))
    for k, e in zip(SWEEP_K, SWEEP_EX):
        print('%12.4f %12.3f' % (k, e))
    print('  E_z was 7.493 / 7.487 / 7.476 -- unchanged, as a layer parallel')
    print('  to the load is indifferent to the fill\'s bulk modulus.')

    k_need = K_for(E_MEASURED)
    print('\nK that reproduces the measured base modulus %.2f GPa : %.4f GPa'
          % (E_MEASURED, k_need))
    print('that is %.1f%% of pure brine.' % (100.0 * k_need / K_BRINE))

    print('\nWood\'s law: free gas needed in the brine to reach that K')
    print('%12s %14s %12s' % ('gas in fill', 'K_mix GPa', 'E_x GPa'))
    for f in (0.0, 1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 4e-2):
        km = float(wood(f))
        # only interpolate inside the sweep, do not extrapolate a modulus
        ex = ('%12.3f' % np.exp(np.interp(np.log(km),
              np.log(SWEEP_K[::-1]), np.log(SWEEP_EX[::-1])))
              if SWEEP_K.min() <= km <= SWEEP_K.max() else '%12s' % '--')
        print('%12.4f %14.4f %s' % (f, km, ex))

    # invert Wood for the gas fraction matching k_need
    fs = np.logspace(-6, -1, 20001)
    f_need = float(fs[int(np.argmin(np.abs(wood(fs) - k_need)))])
    print('\ngas fraction of the fill that lands on the measured modulus:'
          ' %.3f%%' % (100.0 * f_need))

    print('\nFor scale: the decks carry a gas fraction of the CELL of about')
    print('0.006-0.020, against a brine fraction of 0.15-0.23. If that gas sat')
    print('inside the brine rather than in separate voids, the fill would hold')
    gc, bc = 0.010, 0.19
    f_real = gc / (gc + bc)
    print('f = %.3f of free gas, giving K_mix = %.5f GPa -- well BELOW what is'
          % (f_real, wood(f_real)))
    print('needed. So the mechanism has ample room; the modulus is set by how')
    print('much of the gas is co-located with the brine, which is a')
    print('microstructural question and a measurable one.')

    print('\nWHAT THIS DOES NOT SETTLE. Co-locating gas is one way to soften')
    print('the fill; drainage is another, and a slow test on permeable ice')
    print('relaxes the pore pressure to the same effect. Both give a compliant')
    print('fill and this calculation cannot tell them apart -- it only shows')
    print('the required compliance is physically ordinary rather than extreme.')


if __name__ == '__main__':
    main()
