# -*- coding: utf-8 -*-
"""Is the soft base Kujala reports a property of the ice, or of their inversion?

Kujala et al. (1990) do not measure E(z). They measure two things per beam --
the deflection under load, which fixes the flexural rigidity, and the strains
at the two surfaces, whose ratio fixes the neutral-axis position -- and then
solve for the two parameters of a profile ASSUMED LINEAR in depth. That is a
sound reduction if E(z) really is linear. Ours is not: it is convex, falling
slowly through the cold upper column and steeply near the base.

Fitting a straight line to the observables of a convex profile does not return
the convex profile's endpoints. The question is by how much it misses, and in
which direction. If a linear inversion of OUR OWN profile returns a bottom
modulus far below our actual bottom modulus, then a good part of the
disagreement is a like-for-like failure in the comparison rather than missing
knockdown in the model, and chasing the gap with more morphology would be
chasing an artefact.

The test inverts our profile exactly the way they inverted theirs:

    D  = int E(z) (z - z0)^2 dz      flexural rigidity, from deflection
    z0 : int E(z) (z - z0) dz = 0    neutral axis, from the strain ratio

Two observables, two unknowns, so a unique linear profile reproduces both. It
is then compared against the true profile it was computed from.
"""
import csv
import os
import re
import sys

import numpy as np

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = np.array([0.37, 0.38, 0.39, 0.38])


def neutral_axis(E, z):
    """z0 solving int E (z - z0) dz = 0."""
    return np.trapz(E * z, z) / np.trapz(E, z)


def rigidity(E, z, z0):
    return np.trapz(E * (z - z0) ** 2, z)


def linear_profile(Et, Eb, z, H):
    return Et + (Eb - Et) * (z / H)


def invert_linear(D_target, z0_target, z, H):
    """Find (Et, Eb) whose linear profile has this rigidity and neutral axis.

    z0 depends only on the SHAPE, so on the ratio r = Eb/Et; solve that first
    by bisection, then scale to match D. Exactly the two-observable reduction
    the field measurement performs.
    """
    def z0_of_ratio(r):
        E = linear_profile(1.0, r, z, H)
        return neutral_axis(E, z)

    lo, hi = 1e-4, 1.0
    # z0/H rises with r; bracket then bisect
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if z0_of_ratio(mid) < z0_target:
            lo = mid
        else:
            hi = mid
    r = 0.5 * (lo + hi)
    E = linear_profile(1.0, r, z, H)
    z0 = neutral_axis(E, z)
    scale = D_target / rigidity(E, z, z0)
    return scale, scale * r, r


def load_profile(path, pcol, ecol):
    d, e = [], []
    for row in csv.DictReader(open(path, encoding='utf8', errors='replace')):
        m = re.search(r'z(\d+)', row.get('run_id', ''))
        if not m:
            continue
        try:
            v = float(row[ecol])
        except (ValueError, TypeError, KeyError):
            continue
        d.append(float(m.group(1)))
        e.append(v / 1e9)
    if not d:
        return None, None
    d = np.array(d, dtype=float)
    d = d / d.max()
    e = np.array(e)
    zs = sorted(set(d))
    return (np.array(zs),
            np.array([e[d == z].mean() for z in zs]))


def main():
    root = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), 'results')
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        root, 'results_gogo_column.csv')
    z, E = load_profile(src, 'phi_soft_total', 'E_x')
    if z is None:
        print('no usable profile in %s' % src)
        return 1

    H = 1.0
    z = z * H
    z0 = neutral_axis(E, z)
    D = rigidity(E, z, z0)

    print('TRUE computed profile (%s)' % os.path.basename(src))
    print('  E(top)  = %.3f GPa' % E[0])
    print('  E(base) = %.3f GPa' % E[-1])
    print('  z0/H    = %.3f' % (z0 / H))
    print('  D       = %.4f GPa m^3' % D)

    Et, Eb, r = invert_linear(D, z0, z, H)
    print('\nSame beam, reduced the way the field data are reduced')
    print('(linear E(z) fitted to rigidity and neutral axis):')
    print('  E(top)  fitted = %.3f GPa   (true %.3f, %+.1f%%)'
          % (Et, E[0], 100 * (Et / E[0] - 1)))
    print('  E(base) fitted = %.3f GPa   (true %.3f, %+.1f%%)'
          % (Eb, E[-1], 100 * (Eb / E[-1] - 1)))
    print('  ratio  Eb/Et   = %.3f' % r)

    print('\nAgainst Kujala\'s four beams')
    print('  their E(top)  %.2f-%.2f GPa  (mean %.2f)'
          % (K_TOP.min(), K_TOP.max(), K_TOP.mean()))
    print('  their E(base) %.2f-%.2f GPa  (mean %.2f)'
          % (K_BOT.min(), K_BOT.max(), K_BOT.mean()))
    print('  their z0/H    %.2f-%.2f' % (K_Z0.min(), K_Z0.max()))

    print('\nWHAT THIS DOES AND DOES NOT SHOW')
    gap_raw = E[-1] / K_BOT.mean()
    gap_fit = Eb / K_BOT.mean()
    print('  base modulus, true profile vs their mean : %.2fx too stiff' % gap_raw)
    print('  base modulus, inverted the same way      : %.2fx too stiff' % gap_fit)
    if gap_fit < gap_raw:
        print('  so the inversion accounts for a factor %.2f of the gap,'
              % (gap_raw / gap_fit))
        print('  leaving %.2fx that is genuinely missing knockdown.' % gap_fit)
    else:
        print('  the inversion does NOT explain the gap; it is real.')
    print('\n  The neutral axis is the sharper test, being independent of any')
    print('  overall stiffness scale: ours inverts to z0/H = %.3f against their'
          % (neutral_axis(linear_profile(Et, Eb, z, H), z) / H))
    print('  measured %.2f-%.2f. A profile that cannot reproduce THAT is the'
          % (K_Z0.min(), K_Z0.max()))
    print('  wrong shape, whatever its endpoints.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
