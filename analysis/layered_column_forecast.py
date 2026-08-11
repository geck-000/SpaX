# -*- coding: utf-8 -*-
r"""If the basal knockdown is layered, what does the whole column then look like?

The layered cells give a transverse modulus at real plate spacing. The question
this answers is whether wiring that into the column moves the three quantities
the field data actually constrain -- the base modulus, the grading parameter
alpha = E_base/E_top, and the neutral-axis position -- toward the measurements
or past them.

Morphology is not uniform down the column and must not be applied as if it
were. Cold ice holds brine in isolated pockets and warm ice in connected
sheets, which is an observation \citep{light_effects_2003} and not a free
choice, so the layered law is switched on with depth. The switch is tied to the
percolation threshold, which is where connectivity appears, rather than fitted.

Everything here is an ESTIMATE pending the layered sweep across phi at physical
spacing. The layered moduli are known at phi = 0.15 for several spacings and
across phi only at one spacing, so the two are combined multiplicatively. That
is good enough to answer 'closer or further', which is the question, and not
good enough to quote.
"""
import numpy as np

E_ICE = 9.37

# measured targets
K_TOP, K_BOT = 8.05, 1.27          # Kujala, four strain-gauged beams
K_ALPHA = K_BOT / K_TOP
K_Z0 = (0.37, 0.39)
M_ALPHA = 0.38                     # Marchenko / Kerr-Palmer fit

# layered cells, drained, one layer, across phi   (local bracket run)
PHI_L = np.array([0.10, 0.15, 0.227])
E_L_N1 = np.array([1.150, 1.039, 0.724])
# spacing factor from one layer to four, measured at phi = 0.15 on the cluster
SPACING_FACTOR = 0.357 / 0.942
E_U_N1 = np.array([4.705, 2.866, 2.127])
SPACING_FACTOR_U = 5.442 / 2.640


def pocket(phi):
    return E_ICE * (1.0 - 1.65 * phi)


def layered(phi, drained=True):
    """Layered transverse modulus at physical spacing (four layers)."""
    if drained:
        return np.interp(phi, PHI_L, E_L_N1) * SPACING_FACTOR
    return np.minimum(np.interp(phi, PHI_L, E_U_N1) * SPACING_FACTOR_U, E_ICE)


def column(w, phi, z, drained=True):
    """Blend the two morphologies in log-modulus, w = layered share."""
    Ep = pocket(phi)
    El = layered(phi, drained)
    return np.exp((1.0 - w) * np.log(Ep) + w * np.log(np.maximum(El, 1e-6)))


def neutral(E, z):
    return np.trapz(E * z, z) / np.trapz(E, z)


def report(name, E, z):
    a = E[-1] / E[0]
    z0 = neutral(E, z)
    print('%-34s %7.3f %7.3f %7.3f %7.3f'
          % (name, E[0], E[-1], a, z0))
    return a, z0


def main():
    z = np.linspace(0, 1, 400)
    # the column's own C-shaped porosity
    phi = np.interp(z, [0, .29, .63, .79, .96, 1.0],
                    [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])

    print('%-34s %7s %7s %7s %7s'
          % ('', 'E_top', 'E_base', 'alpha', 'z0/H'))
    print('-' * 66)
    report('pockets throughout (current)', pocket(phi), z)

    # layered share switched on at the percolation threshold, over a transition
    for width in (0.25,):
        for phic in (0.05,):
            w = 1.0 / (1.0 + np.exp(-(phi - phic) / (width * phic)))
            w = (w - w.min()) / (w.max() - w.min())
            report('pocket->layer, drained', column(w, phi, z, True), z)
            report('pocket->layer, undrained', column(w, phi, z, False), z)

    # layers only in the warm base, a resolved skeletal zone
    for zc in (0.85, 0.75):
        w = np.clip((z - zc) / (1.0 - zc), 0.0, 1.0)
        report('layers below z/H = %.2f, drained' % zc,
               column(w, phi, z, True), z)

    print('-' * 66)
    print('%-34s %7.2f %7.2f %7.3f %7s'
          % ('Kujala (measured)', K_TOP, K_BOT, K_ALPHA,
             '%.2f-%.2f' % K_Z0))
    print('%-34s %7s %7s %7.3f %7s'
          % ('Marchenko / Kerr-Palmer', '-', '-', M_ALPHA, '-'))

    print('\nREADING THIS')
    print('  alpha is the sharper of the two, being a ratio: it does not care')
    print('  about any overall stiffness offset. z0/H is sharper still, since')
    print('  it is fixed by the SHAPE alone. A model that gets the base modulus')
    print('  by softening only the last few percent of the thickness will move')
    print('  alpha and leave z0 where it was, which is the failure the current')
    print('  column already has.')
    print('\n  Estimates only: the layered law is known across phi at one')
    print('  spacing and across spacing at one phi, and is combined')
    print('  multiplicatively here. The campaign that removes that assumption')
    print('  is rve_bracket_layer at physical spacing.')


if __name__ == '__main__':
    main()
