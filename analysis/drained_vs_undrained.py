# -*- coding: utf-8 -*-
"""How much of the column's stiffness is the incompressible-brine assumption?

The cells model brine as a soft solid with K = 2.2 GPa and G = 4.4e-4 GPa: a
liquid, carrying no shear but resisting compression like water. That is the
UNDRAINED response, correct for loading fast enough that brine cannot leave the
pocket it is in.

The measurements it is being compared against are not that. Kujala's four-point
beams and Gogolaze's cantilevers are loaded over seconds to minutes, on warm
permeable ice connected to the ocean, where brine can move. The relevant
response there is drained: the pocket pressure relaxes and the brine carries
almost nothing.

Mori-Tanaka is used rather than the cells because the question is about the
phase moduli and not about geometry, and MT is exact to first order in phi for
either filling. If the drained-undrained gap is small the incompressibility
assumption is not the problem and the search goes back to topology; if it is
large it has to be settled before any topology conclusion means anything.
"""
import numpy as np

E_ICE, NU_ICE = 9.37, 0.33
K_M = E_ICE / (3.0 * (1.0 - 2.0 * NU_ICE))
G_M = E_ICE / (2.0 * (1.0 + NU_ICE))

# brine as modelled (undrained) and as a drained skeleton
K_UNDRAINED, G_BRINE = 2.2, 4.4e-4
K_DRAINED = 2.2e-3


def mori_tanaka(phi, K_i, G_i, K_m=K_M, G_m=G_M):
    """Spherical inclusions at volume fraction phi."""
    alpha = 3.0 * K_m / (3.0 * K_m + 4.0 * G_m)
    beta = 6.0 * (K_m + 2.0 * G_m) / (5.0 * (3.0 * K_m + 4.0 * G_m))
    # dilute strain-concentration factors for a sphere
    a = 1.0 / (1.0 + alpha * (K_i / K_m - 1.0))
    b = 1.0 / (1.0 + beta * (G_i / G_m - 1.0))
    K = K_m + phi * (K_i - K_m) * a / (1.0 - phi + phi * a)
    G = G_m + phi * (G_i - G_m) * b / (1.0 - phi + phi * b)
    E = 9.0 * K * G / (3.0 * K + G)
    return K, G, E


def main():
    print('ice: K = %.3f, G = %.3f GPa\n' % (K_M, G_M))
    print('%8s %11s %11s %9s   %s' % (
        'phi', 'E undrained', 'E drained', 'ratio', 'note'))
    for phi in (0.05, 0.086, 0.104, 0.15, 0.20, 0.227, 0.30):
        _, _, eu = mori_tanaka(phi, K_UNDRAINED, G_BRINE)
        _, _, ed = mori_tanaka(phi, K_DRAINED, G_BRINE)
        note = ''
        if abs(phi - 0.086) < 1e-9:
            note = 'column minimum'
        elif abs(phi - 0.104) < 1e-9:
            note = 'column top'
        elif abs(phi - 0.227) < 1e-9:
            note = 'column base'
        print('%8.3f %11.3f %11.3f %9.3f   %s' % (phi, eu, ed, eu / ed, note))

    print('\nwhat the column would do, top to base')
    _, _, u_top = mori_tanaka(0.104, K_UNDRAINED, G_BRINE)
    _, _, u_bot = mori_tanaka(0.227, K_UNDRAINED, G_BRINE)
    _, _, d_top = mori_tanaka(0.104, K_DRAINED, G_BRINE)
    _, _, d_bot = mori_tanaka(0.227, K_DRAINED, G_BRINE)
    print('  undrained  %.3f -> %.3f GPa   contrast %.2f' % (
        u_top, u_bot, u_top / u_bot))
    print('  drained    %.3f -> %.3f GPa   contrast %.2f' % (
        d_top, d_bot, d_top / d_bot))
    print('  measured   %.2f -> %.2f GPa   contrast %.2f' % (
        8.05, 1.27, 8.05 / 1.27))

    print('\nA column that is undrained at the top and drained at the base --')
    print('the permeability transition -- would span')
    print('  %.3f -> %.3f GPa   contrast %.2f' % (u_top, d_bot, u_top / d_bot))
    print('against the measured %.2f. Drainage alone moves the base by a factor'
          % (8.05 / 1.27))
    print('%.2f, so it is a real term but not the whole gap.' % (u_bot / d_bot))


if __name__ == '__main__':
    main()
