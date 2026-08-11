r"""How close are we to the stated goal: E(z) at a reasonable level and shape?

The goal is not to reproduce the strength-style correlations. Those sit below
the Hashin-Shtrikman lower bound -- Weeks and Assur give 2.37 GPa at phi = 0.086
against a bound of 7.99 -- so they contain something a linear elastic two-phase
model does not, most plausibly damage. Chasing them would mean fitting damage
into an elastic closure. The target instead is the depth profile: a modulus
that starts near the measured cold-end value, falls monotonically, and has
roughly the right grading and curvature.

Three things define that target and each is checked separately:

  level     E at the cold end, where Kujala's four beams give 7.18-8.60 GPa
  grading   alpha = E_base / E_top; Marchenko's Kerr-Palmer fit gives 0.38
  shape     curvature: measured profiles fall fastest at the surface and
            flatten with depth, so the profile sits BELOW its own chord

The porosity profile is treated as an input to be stated, not as part of the
microstructure model, because it is where our synthetic column and the measured
ice differ most.
"""
import numpy as np

import layered_law as law
from shape_diagnosis import marchenko_E, corr_inv, ours_phi

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
M_ALPHA, M_Z0 = 0.384, 0.428


def metrics(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    En = E / E[0]
    chord = En[0] + (En[-1] - En[0]) * z
    return E[0], En[-1], z0, float(np.mean(En - chord))


def main():
    z = np.linspace(1e-3, 1.0, 400)
    his = corr_inv(marchenko_E(z))          # monotonic, 0.014 -> 0.122
    our = ours_phi(z)                       # C-shaped, 0.104 -> 0.086 -> 0.227

    print('TARGET  E_top 7.18-8.60 GPa (Kujala) | alpha 0.38 | z0/H 0.43 |'
          ' curvature negative')
    print()
    for tag, phi in (('a MONOTONIC phi (Marchenko implied)', his),
                     ('OUR synthetic C-shaped phi', our)):
        print('--- on %s ---' % tag)
        print('%-26s %8s %8s %8s %10s' % ('', 'E_top', 'alpha', 'z0/H', 'curv'))
        for name, E in (('pockets only', law.pocket(phi)),
                        ('layered, b^1', law.layered(phi, 1.0)),
                        ('layered, b^2', law.layered(phi, 2.0))):
            t, a, z0, c = metrics(E, z)
            flag = ''
            if K_TOP.min() <= t <= K_TOP.max():
                flag += ' level ok'
            if abs(a - M_ALPHA) < 0.12:
                flag += ' grading ok'
            if c < 0:
                flag += ' shape ok'
            print('%-26s %8.2f %8.3f %8.3f %+10.3f  %s'
                  % (name, t, a, z0, c, flag))
        print()

    print('READING')
    print('  On a monotonic porosity the layered closure meets all three at')
    print('  once, and the exponent barely matters for the SHAPE -- both b^1')
    print('  and b^2 give negative curvature. It matters for the grading.')
    print('  On our C-shaped porosity no exponent gives negative curvature,')
    print('  because a modulus that follows a C-shaped porosity must rise')
    print('  somewhere. The obstacle to the stated goal is the porosity')
    print('  profile, not the microstructure model.')
    print()
    print('  That is a good position to be in: phi(z) is an INPUT, set by the')
    print('  temperature and salinity profile a user supplies, and it is')
    print('  reported rather than fitted. The microstructure model should be')
    print('  validated in phi, where it is being asked a question it can')
    print('  answer, and the column then follows from whatever phi(z) the ice')
    print('  in question actually has.')


if __name__ == '__main__':
    main()
