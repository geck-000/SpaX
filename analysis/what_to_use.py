r"""What the evidence actually supports, and what each option costs.

MEASURED BY OUR OWN CELLS, and not in dispute:

  pocket knockdown        E = 9.37 (1 - 1.65 phi), R^2 = 0.999, and Mori-Tanaka
                          reproduces it to under 1% with nothing fitted
  drainage, pockets       1.04x between the limits: irrelevant
  drainage, layers        2.8-6.9x: decisive, because a confined layer resists
                          at its bulk modulus
  undrained layers        a genuine RVE, CV 0.9% at fixed spacing
  bridge subdivision      E ~ N^0.458 at fixed area, i.e. constriction (0.500),
                          not bending and not area
  plate spacing           E ~ a0^0.69 drained, at fixed phi, b and bridge count

MEASURED BY OTHERS, and usable:

  a0 = 200-500 um         Pringle, direct from tomography
  phi(T)                  Pringle Table 1, and it validates Frankenstein-Garner
                          to +3.2% once his CsCl doping is corrected for

ASSUMED, AND NOT SECURE:

  b = 1 - sqrt(phi/phi_0) Assur's geometry -- which Pringle's own paper says is
                          "much more complicated than suggested by simple
                          models of parallel ice lamellae"
  phi_0 = 0.20            defensible range 0.15-0.36, and the modulus is
                          roughly linear in b
  nu                      not obtainable: published inclusion densities differ
                          24x through imaging resolution alone

So the layered closure has one measured exponent, one measured length, and two
parameters that cannot be pinned. That is too many to call it predictive, and
too few to call it useless.
"""
import numpy as np

import layered_law as law

EXPONENTS = {'our cells (constriction)': 0.50,
             'Gogolaze at a0 = 0.35 mm': 0.62,
             'Marchenko (spacing-independent)': 0.83}
A0_MM, A0_REF_MM, SPACING_EXP = 0.35, 0.75, 0.69


def layered(phi, n):
    b = np.maximum(law.assur_b(phi), 0.0)
    return law.pocket(phi) * b ** n * (A0_MM / A0_REF_MM) ** SPACING_EXP


def main():
    ns = np.array(list(EXPONENTS.values()))
    print('THE THREE EXPONENTS, once the measured spacing is used')
    for k, v in EXPONENTS.items():
        print('  %-34s %.2f' % (k, v))
    print('  spread %.2f to %.2f, mid %.2f' % (ns.min(), ns.max(), ns.mean()))

    print('\nWHAT THAT SPREAD COSTS IN MODULUS')
    print('%8s %10s %10s %10s %10s' % ('phi', 'b', 'n=0.50', 'n=0.83', 'spread'))
    for phi in (0.05, 0.10, 0.15, 0.19):
        b = float(law.assur_b(phi))
        lo, hi = float(layered(phi, 0.83)), float(layered(phi, 0.50))
        print('%8.3f %10.3f %10.3f %10.3f %9.0f%%'
              % (phi, b, hi, lo, 100 * (hi / lo - 1)))

    print('\n  So the exponent uncertainty alone is worth 30-90% in the layered')
    print('  modulus, before phi_0 and nu are counted. phi_0 over its own range')
    print('  of 0.15-0.36 moves b by more again.')

    print('\nTHE OPTIONS')
    print('  A. Report the pocket column as the quantitative result, and the')
    print('     layered work as a mechanism that BOUNDS the base. Everything')
    print('     quoted is then measured. The base stays a bound, which is what')
    print('     the paper already says, but now with the mechanism named and')
    print('     four of its ingredients measured rather than asserted.')
    print('\n  B. Quote the layered closure with n = %.2f and state the'
          % ns.mean())
    print('     uncertainty honestly. Gets a usable E(z) at the cost of two')
    print('     parameters that cannot currently be measured.')
    print('\n  C. Keep pockets above the porosity where the plane empties and')
    print('     an explicitly EMPIRICAL knockdown below it. Honest, usable, and')
    print('     makes no claim the microstructure supports.')

    print('\nRECOMMENDATION: A, with B in a clearly separated subsection.')
    print('  The pocket result is verified and nothing in this work has')
    print('  weakened it. The layered result is a real mechanism -- confinement,')
    print('  constriction and spacing are all measured, and two of them were')
    print('  measured against predictions that could have failed and did not.')
    print('  What it is not yet is a closure, because b and nu are not pinned,')
    print('  and Pringle says the geometry b comes from is too simple for his')
    print('  images. Presenting it as a mechanism is defensible; presenting it')
    print('  as a predictive E(z) would need those two measurements first.')


if __name__ == '__main__':
    main()
