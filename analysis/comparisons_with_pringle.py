r"""The comparisons redone with Pringle's measured plate spacing.

Two things from the paper enter, and only two, because they are the only two it
actually supplies for our purposes.

  a0 = 200-500 um, the ice lamella thickness, measured directly. Every campaign
  we ran used 0.75-3 mm, so the physical spacing is finer than anything solved.

  phi(T) at S = 9.4 ppt, which after correcting for the CsCl doping validates
  Frankenstein-Garner to +3.2% and so needs no change to the column.

Spacing matters because the constrictions add in series down the column: a
layer plane per a0 of thickness means 1/a0 constrictions per unit length, so
finer lamellae are SOFTER. Our own layer-count sweep measured that directly at
fixed porosity, bridge fraction and bridge density, giving E ~ a0^0.69 for the
drained cell.

So the closure carries a spacing factor,

    E_layer(phi) = E_pocket(phi) * b(phi)^n * (a0 / a0_ref)^0.69

with a0_ref the spacing the exponent n was fitted at. Putting the measured a0
in and refitting n is the test: if the exponents move onto the constriction
value the cells themselves report, three independent routes agree.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

import layered_law as law

SPACING_EXP = 0.69           # measured, drained, layer-count sweep
CELL_MM = 3.0                # our cell edge
A0_REF_MM = CELL_MM / 4.0    # n_slabs = 4, the finest we solved
A0_PRINGLE = (0.20, 0.50)    # measured ice lamella thickness, mm

GOGO_COR, H_GOGO = 1.421, 0.32
M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6
C_A, C_B = 7.23, 4.2


def marchenko_E(d):
    return M_E0 * (1.0 - (1.0 - M_ALPHA) * d ** M_N)


def corr_inv(E):
    return (np.log(C_A / np.maximum(E, 1e-9)) / C_B) ** 2


def layered(phi, n, a0_mm):
    return law.pocket(phi) * law.assur_b(phi) ** n * \
        (a0_mm / A0_REF_MM) ** SPACING_EXP


def flexural(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def main():
    z = np.linspace(1e-3, 1.0, 400)
    print('SPACING FACTOR from Pringle, relative to the finest cells we solved')
    print('  a0_ref = %.2f mm (n_slabs = 4 in a %.0f mm cell)'
          % (A0_REF_MM, CELL_MM))
    for a0 in A0_PRINGLE:
        print('  a0 = %.2f mm -> factor %.3f, i.e. %.0f%% softer'
              % (a0, (a0 / A0_REF_MM) ** SPACING_EXP,
                 100 * (1 - (a0 / A0_REF_MM) ** SPACING_EXP)))

    zc = z * H_GOGO * 100.0
    phi_g = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    phi_m = corr_inv(marchenko_E(z))
    tgt_m = marchenko_E(z)

    print('\nEXPONENT EACH DATASET NEEDS, at each measured spacing')
    print('%-16s %10s %10s %10s' % ('a0 (mm)', 'Gogolaze', 'Marchenko', 'note'))
    for a0 in (A0_REF_MM, A0_PRINGLE[1], 0.35, A0_PRINGLE[0]):
        try:
            ng = brentq(lambda n: flexural(layered(phi_g, n, a0), z) - GOGO_COR,
                        0.02, 12.0)
        except ValueError:
            ng = float('nan')

        def miss(n):
            E = layered(phi_m, n, a0)
            return float(np.mean((E / E[0] - tgt_m / tgt_m[0]) ** 2))
        nm = minimize_scalar(miss, bounds=(0.05, 8.0), method='bounded').x
        tag = 'our cells' if abs(a0 - A0_REF_MM) < 1e-9 else 'Pringle range'
        print('%-16.2f %10.2f %10.2f %10s' % (a0, ng, nm, tag))

    print('\n  Marchenko\'s exponent does not move, because his target is a')
    print('  normalised SHAPE and the spacing factor is a constant multiplier')
    print('  that divides out. Gogolaze\'s does move, because his is a level.')

    a0_mid = 0.35
    ng = brentq(lambda n: flexural(layered(phi_g, n, a0_mid), z) - GOGO_COR,
                0.02, 12.0)
    print('\nAT THE MIDDLE OF PRINGLE\'S RANGE, a0 = %.2f mm:' % a0_mid)
    print('  Gogolaze needs b^%.2f' % ng)
    print('  our own cells report constriction, which is b^0.50')
    print('  Marchenko needs b^%.2f'
          % minimize_scalar(lambda n: float(np.mean(
              (layered(phi_m, n, a0_mid) / layered(phi_m, n, a0_mid)[0]
               - tgt_m / tgt_m[0]) ** 2)),
              bounds=(0.05, 8.0), method='bounded').x)

    print('\n  Whether that counts as agreement depends on how close is close,')
    print('  and the spacing factor was extrapolated a factor of two below the')
    print('  finest cell actually solved, on a power law fitted over four')
    print('  points. It is suggestive, not settled.')


if __name__ == '__main__':
    main()
