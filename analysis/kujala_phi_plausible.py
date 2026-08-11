r"""Is the porosity Kujala's beams imply reachable by real ice?

Feeding the implied porosity back into the closure is circular: it was obtained
by inverting his profile THROUGH that closure, so it reproduces his profile
exactly and tests nothing. The question with content is whether a physically
plausible temperature and salinity profile can produce it.

Brine volume follows Frankenstein and Garner,

    phi = S (-49.185 / T + 0.532) / 1000,   T in Celsius, S in ppt,

so for an assumed thermal profile the implied porosity dictates the salinity
profile required. Sea ice bulk salinity is bounded: first-year ice runs roughly
4 to 10 ppt, C-shaped, and cannot exceed the ocean it froze from. If the
implied porosity demands salinities outside that, it is not ice, and the
agreement obtained by using it would be an artefact of the inversion.

Three thermal profiles are tried, since Kujala reports none: a cold Baltic
winter, a mild one, and a near-isothermal spring column.
"""
import numpy as np
from scipy.optimize import brentq

import layered_law as law

K_TOP, K_BOT = 8.05, 1.27
S_OCEANIC = (2.0, 12.0)        # ppt, oceanic first-year sea ice bulk salinity
# Kujala's beams are BALTIC, and the Baltic is brackish: ice grown from it
# carries bulk salinities far below oceanic first-year values. Which bound
# applies changes the verdict entirely, so both are reported.
S_BALTIC = (0.2, 2.0)
S_PLAUSIBLE = S_OCEANIC
PROFILES = (('cold winter  -20 to -1.8 C', -20.0, -1.8),
            ('mild winter  -12 to -1.8 C', -12.0, -1.8),
            ('spring        -5 to -1.8 C', -5.0, -1.8))


def implied_phi(z, exponent):
    tgt = K_TOP + (K_BOT - K_TOP) * z
    out = []
    for t in tgt:
        try:
            out.append(brentq(lambda p: law.layered(p, exponent) - t,
                              1e-6, 0.95))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def salinity_needed(phi, T):
    """Invert Frankenstein-Garner for S at each depth."""
    return phi * 1000.0 / (-49.185 / T + 0.532)


def main():
    z = np.linspace(1e-3, 1.0, 200)
    print('The implied porosity is what it is; the test is what ice it needs.')
    print('First-year bulk salinity runs about %.0f-%.0f ppt.\n'
          % S_PLAUSIBLE)

    for exponent in (0.63, 0.99):
        phi = implied_phi(z, exponent)
        print('=' * 68)
        print('IMPLIED BY KUJALA WITH b^%.2f : phi %.3f -> %.3f'
              % (exponent, phi[0], phi[-1]))
        print('=' * 68)
        print('%-30s %9s %9s %9s %s'
              % ('thermal profile', 'S(top)', 'S(mid)', 'S(base)', 'verdict'))
        for name, Ts, Tb in PROFILES:
            T = Ts + (Tb - Ts) * z
            S = salinity_needed(phi, T)
            smid = float(np.interp(0.5, z, S))
            ok = (S_PLAUSIBLE[0] <= np.nanmin(S) and
                  np.nanmax(S) <= S_PLAUSIBLE[1])
            worst = np.nanmax(S)
            verdict = 'plausible' if ok else 'needs %.0f ppt' % worst
            print('%-30s %9.1f %9.1f %9.1f %s'
                  % (name, S[0], smid, S[-1], verdict))
        print()

    print('=' * 68)
    print('WHICH SALINITY BOUND APPLIES')
    print('=' * 68)
    print('  Oceanic first-year ice  %.1f-%.1f ppt' % S_OCEANIC)
    print('  BALTIC ice              %.1f-%.1f ppt -- and Kujala\'s beams are'
          % S_BALTIC)
    print('  Baltic, grown from brackish water. He reports no salinity, but the')
    print('  basin does not supply oceanic values.')
    for exponent in (0.63, 0.99):
        phi = implied_phi(z, exponent)
        best = min(np.nanmax(salinity_needed(phi, Ts + (Tb - Ts) * z))
                   for _, Ts, Tb in PROFILES)
        print('\n  b^%.2f needs at least %.1f ppt somewhere in the column:'
              % (exponent, best))
        print('     against oceanic %.0f-%.0f ppt : %s'
              % (S_OCEANIC[0], S_OCEANIC[1],
                 'possible' if best <= S_OCEANIC[1] else 'NOT possible'))
        print('     against Baltic  %.1f-%.1f ppt : %s'
              % (S_BALTIC[0], S_BALTIC[1],
                 'possible' if best <= S_BALTIC[1] else
                 'NOT possible, %.0fx over' % (best / S_BALTIC[1])))

    print('\nREADING')
    print('  Frankenstein-Garner has phi proportional to S/|T|, so a porosity')
    print('  that rises steeply with depth can be produced either by salinity')
    print('  rising or by the ice warming. Warming does most of the work near')
    print('  the base, where |T| falls toward 1.8, but not in mid-column, and')
    print('  that is where the required salinity is worst.')
    print('\n  If no profile in the table is plausible, then the porosity that')
    print('  reproduces Kujala under this closure is not a porosity real ice')
    print('  has, and the closure is compensating with an unphysical input')
    print('  rather than describing his beams.')


if __name__ == '__main__':
    main()
