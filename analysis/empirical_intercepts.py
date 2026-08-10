# -*- coding: utf-8 -*-
"""Separating what the empirical E(v_b) laws say about ice from what they say
about how ice was measured.

The laws disagree at zero brine by a factor of three -- Karulina extrapolates to
3.10 GPa, Vaudrey to 5.31, Marchenko to 7.23, Weeks & Assur to 9.5. Pure ice is
about 9.3. A correlation that returns 3 GPa for brine-free ice is not describing
microstructure at that end; it carries the compliance, scale and rate of the
test it was fitted from. Comparing a homogenisation, which is anchored at the
pure-ice value by construction, against the raw curves therefore charges the
model for an offset it should not be asked to reproduce.

Normalising each law by its own intercept removes that, and what is left is the
shape: how fast the modulus falls per unit brine. That residual is the real
target for a microstructural model, and it is where our cells actually fail.

Weeks & Assur is the informative case. It is the only law that extrapolates to
the right pure-ice modulus, so its steepness cannot be dismissed as an
intercept artefact, and its form (1 - sqrt(v))^4 is an Assur load-bearing area
raised to a power near the percolation exponent of a 3D elastic network.
"""
import numpy as np

E_ICE = 9.37


def weeks(v):
    return 9.5 * (1.0 - np.sqrt(v)) ** 4


def vaudrey(v):
    return 5.31 - 0.436 * np.sqrt(v * 1000.0)


def karulina(v):
    return 3.1031 * np.exp(-3.385 * np.sqrt(v))


def marchenko(v):
    return 7.23 * np.exp(-4.2 * np.sqrt(v))


def ours(v):
    return E_ICE * (1.0 - 1.65 * v)


LAWS = [('Weeks & Assur 1967', weeks), ('Vaudrey 1977', vaudrey),
        ('Karulina 2019', karulina), ('Marchenko 2024', marchenko),
        ('SpaX cells', ours)]


def main():
    print('intercept at v_b = 0   (pure ice is %.2f GPa)\n' % E_ICE)
    for name, f in LAWS:
        e0 = float(f(1e-12))
        print('  %-20s %6.2f GPa   %+.0f%% vs ice' % (
            name, e0, 100.0 * (e0 / E_ICE - 1.0)))

    print('\nraw values, GPa')
    vs = (0.02, 0.05, 0.09, 0.15, 0.227)
    print('  %-20s' % 'v_b' + ''.join('%9.3f' % v for v in vs))
    for name, f in LAWS:
        print('  %-20s' % name + ''.join('%9.3f' % f(v) for v in vs))

    print('\nnormalised by each law\'s own intercept -- shape only')
    print('  %-20s' % 'v_b' + ''.join('%9.3f' % v for v in vs))
    for name, f in LAWS:
        e0 = float(f(1e-12))
        print('  %-20s' % name + ''.join('%9.3f' % (f(v) / e0) for v in vs))

    print('\nhow much steeper than our cells, after both are normalised')
    print('  %-20s' % 'v_b' + ''.join('%9.3f' % v for v in vs))
    o0 = ours(1e-12)
    for name, f in LAWS[:-1]:
        e0 = float(f(1e-12))
        print('  %-20s' % name + ''.join(
            '%9.2f' % ((ours(v) / o0) / (f(v) / e0)) for v in vs))

    print('\nWeeks & Assur is the law to answer to: right intercept, so its')
    print('steepness is microstructure and not test compliance. At the column')
    print('base it asks for %.2f GPa where our cells give %.2f.'
          % (weeks(0.227), ours(0.227)))
    print('\nIts form is (1 - sqrt(v))^4. The bracket is Assur load-bearing')
    print('area; the exponent 4 is close to the elastic percolation exponent')
    print('of a 3D network, which is what a solid skeleton being severed by')
    print('brine sheets should follow. Our cells never approach that threshold:')
    print('the matrix stays fully connected at every porosity we build, so they')
    print('sit near the Hashin-Shtrikman upper bound and fall linearly.')

    print('\nCAVEAT, and it cuts both ways. These laws were fitted over field')
    print('brine volumes of roughly 0.05-0.20, so v_b = 0 is an extrapolation')
    print('for all of them and none is evidence about brine-free ice. The')
    print('intercept spread is a reason to distrust the raw curves at the cold')
    print('end, but the normalised shapes inherit that same uncertainty and are')
    print('not clean either. What survives both objections is the comparison')
    print('INSIDE the fitted range: at v_b = 0.05 to 0.23 every law gives')
    print('%.2f-%.2f GPa where our cells give %.2f-%.2f. That gap is not an'
          % (min(f(0.227) for _, f in LAWS[:-1] if f(0.227) > 0),
             max(f(0.05) for _, f in LAWS[:-1]), ours(0.227), ours(0.05)))
    print('artefact of extrapolation, and it is what needs explaining.')
    print('\nA second caveat on the shapes: a 40-45% loss at v_b = 0.02 is')
    print('below the Hashin-Shtrikman lower bound for an isotropic two-phase')
    print('mixture. That is admissible only for a strongly layered or')
    print('connected morphology -- or it means the correlations carry damage')
    print('and cracking that a linear-elastic homogenisation should not chase.')


if __name__ == '__main__':
    main()
