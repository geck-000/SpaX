r"""Why our E(z) has a knee where the measured profiles have a smooth curve.

Two candidate causes, and they call for completely different responses:

  (i)  the morphology switch. Layers turn on at the -5 C isotherm, which puts a
       kink in E(z) that no measured profile shows.
  (ii) the porosity profile. Our column is synthetic -- a designed winter
       first-year profile -- and if its phi(z) has the wrong shape then no E(phi)
       law, smooth or otherwise, can produce the right E(z).

These are separable without any new solves. Marchenko does not measure a
modulus profile: he measures brine content, pushes it through an empirical
correlation (his formula 5), and fits the result with a Kerr-Palmer form (his
eq. 17). Inverting the correlation on the fit therefore recovers the brine
profile his published curve actually stands for, and it can be compared with
ours directly, with no fitted scalar anywhere.

If his phi(z) resembles ours, the knee is our morphology model's fault. If it
does not, the shape disagreement was inherited from the assumed thermal and
salinity profiles and was never about microstructure at all.
"""
import numpy as np

import layered_law as law

M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6      # Marchenko eq. (17)
C_A, C_B = 7.23, 4.2                     # his formula (5): E = A exp(-B sqrt(v))


def marchenko_E(d):
    return M_E0 * (1.0 - (1.0 - M_ALPHA) * d ** M_N)


def corr_inv(E):
    """His correlation, inverted: what brine content does this modulus mean?"""
    return (np.log(C_A / np.maximum(E, 1e-9)) / C_B) ** 2


def ours_phi(z):
    return np.interp(z, [0, .29, .63, .79, .96, 1.0],
                     [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])


def shape(E, z):
    """Normalised profile plus the curvature sign that distinguishes them."""
    En = E / E[0]
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    # a straight line between the endpoints; positive means the profile bulges
    # ABOVE the chord, i.e. falls slowly at first and steeply later
    chord = En[0] + (En[-1] - En[0]) * z
    return En, z0, float(np.mean(En - chord))


def main():
    z = np.linspace(1e-3, 1.0, 400)

    print('THE BRINE PROFILE MARCHENKO\'S CURVE STANDS FOR')
    print('(his Kerr-Palmer fit, inverted through his own correlation)')
    print('%8s %12s %12s %10s' % ('z/H', 'his E (GPa)', 'his phi', 'our phi'))
    for zz in (0.0, 0.25, 0.5, 0.75, 1.0):
        e = marchenko_E(max(zz, 1e-9))
        print('%8.2f %12.3f %12.4f %10.4f'
              % (zz, e, corr_inv(e), ours_phi(zz)))

    hp = corr_inv(marchenko_E(z))
    op = ours_phi(z)
    print('\n  his phi rises monotonically %.4f -> %.4f' % (hp[0], hp[-1]))
    print('  ours is C-shaped %.4f -> %.4f -> %.4f'
          % (op[0], op.min(), op[-1]))
    print('  and ours is %.1fx his at the surface.' % (op[0] / hp[0]))

    print('\n' + '=' * 66)
    print('SO WHOSE FAULT IS THE KNEE? Apply OUR laws to HIS porosity.')
    print('=' * 66)
    print('%-34s %8s %8s %9s' % ('', 'alpha', 'z0/H', 'curvature'))
    for name, E in (
            ('Marchenko, his own curve', marchenko_E(z)),
            ('our pockets, on HIS phi', law.pocket(hp)),
            ('our closure, on HIS phi', law.column(z, hp, exponent=2.0)),
            ('our pockets, on OUR phi', law.pocket(op)),
            ('our closure, on OUR phi', law.column(z, op, exponent=2.0))):
        En, z0, cur = shape(E, z)
        print('%-34s %8.3f %8.3f %+9.3f' % (name, En[-1], z0, cur))

    print('\n  Curvature sign is the discriminator: Marchenko falls fastest at')
    print('  the SURFACE and flattens with depth, so his curvature is negative.')
    print('  The switched closure is strongly POSITIVE even on his porosity,')
    print('  which is the knee showing up as a number.')

    print('\n%-34s %8s %8s %9s' % ('WITHOUT THE MORPHOLOGY SWITCH',
                                   'alpha', 'z0/H', 'curvature'))
    for name, E in (('layers everywhere, on HIS phi', law.layered(hp, 2.0)),
                    ('layers everywhere, on OUR phi', law.layered(op, 2.0))):
        En, z0, cur = shape(E, z)
        print('%-34s %8.3f %8.3f %+9.3f' % (name, En[-1], z0, cur))

    print('\nCONCLUSION -- two separate faults, and the switch is one of them.')
    print('  The knee IS the morphology switch. On his own porosity the')
    print('  switched closure gives curvature +0.198 against his -0.073, while')
    print('  the same closure applied throughout gives -0.033: the right sign,')
    print('  with alpha 0.449 against 0.384 and z0/H 0.440 against 0.428.')
    print('  The switch was never well motivated. Columnar ice carries its')
    print('  lamellar substructure at every depth, set by the growth process;')
    print('  what temperature changes is the brine volume in the lamellae, not')
    print('  whether they exist. And no switch is needed, because b = 1 -')
    print('  sqrt(phi) tends to one as phi falls, so a cold layer stops')
    print('  softening of its own accord.')
    print('\n  The second fault is ours and this does not fix it: the synthetic')
    print('  phi(z) is C-shaped and seven times too rich at the surface, and on')
    print('  it even the unswitched closure gives +0.149. No E(phi) law turns a')
    print('  C-shaped porosity into a monotonic modulus profile.')


if __name__ == '__main__':
    main()
