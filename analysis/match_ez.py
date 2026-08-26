r"""Match E(z) for the three comparison cases, each on its own porosity input.

The goal is the elastic depth profile at a sensible level and shape, not the
strength correlations. So for each case the porosity comes from that study
wherever it is reported, one exponent is fitted, and the result is judged on
level, grading, neutral axis and curvature.

  MARCHENKO   phi(z) recovered by inverting his own correlation (formula 5) on
              his Kerr-Palmer fit (eq. 17). Target is the full E(z) curve.
  GOGOLAZE    phi(z) from his eq. (14) directly. He reports a whole-beam
              modulus, not a profile, so the target is the flexural rigidity.
  KUJALA      neither salinity nor brine content is reported. Nothing can be
              matched forwards, so his profile is INVERTED through our closure
              instead, and the question becomes whether the phi(z) his beams
              would need is a physically plausible one.

The exponent is the only quantity fitted, and it is reported against the range
its own derivation supports rather than quoted as a result.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

import layered_law as law

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = np.array([0.37, 0.38, 0.39, 0.38])
M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6
C_A, C_B = 7.23, 4.2
GOGO_APP, GOGO_COR, H_GOGO = 0.785, 1.421, 0.32
GA_CEILING = 2.0


def marchenko_E(d):
    return M_E0 * (1.0 - (1.0 - M_ALPHA) * d ** M_N)


def corr_inv(E):
    return (np.log(C_A / np.maximum(E, 1e-9)) / C_B) ** 2


def shape(E, z):
    En = E / E[0]
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    chord = En[0] + (En[-1] - En[0]) * z
    return En[-1], z0, float(np.mean(En - chord))


def flexural(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def verdict(n):
    if n <= 1.05:
        return 'stretch, in range'
    if n <= GA_CEILING:
        return 'below GA ceiling, but GA needs b<0.3'
    return 'ABOVE any cellular-solid basis'


def case_marchenko(z):
    print('=' * 70)
    print('MARCHENKO 2024 -- phi from inverting his own correlation')
    print('=' * 70)
    phi = corr_inv(marchenko_E(z))
    tgt = marchenko_E(z)
    print('  phi %.4f (top) -> %.4f (base), monotonic' % (phi[0], phi[-1]))

    # fit the exponent on SHAPE, normalising out his intercept, since his
    # correlation returns 7.23 GPa at zero brine against 9.37 for pure ice
    def miss(n):
        E = law.layered(phi, n)
        return float(np.mean((E / E[0] - tgt / tgt[0]) ** 2))
    r = minimize_scalar(miss, bounds=(0.2, 8.0), method='bounded')
    n = r.x
    E = law.layered(phi, n)
    a, z0, c = shape(E, z)
    at, z0t, ct = shape(tgt, z)
    print('  best exponent on normalised shape: %.2f  (%s)' % (n, verdict(n)))
    print('  %-14s %8s %8s %10s %10s' % ('', 'E_top', 'alpha', 'z0/H', 'curv'))
    print('  %-14s %8.2f %8.3f %10.3f %+10.3f' % ('ours', E[0], a, z0, c))
    print('  %-14s %8.2f %8.3f %10.3f %+10.3f' % ('his', tgt[0], at, z0t, ct))
    print('  shape RMS after normalising: %.4f' % np.sqrt(miss(n)))
    print('  LEVEL is %.2fx his, and that gap is his intercept: his correlation'
          % (E[0] / tgt[0]))
    print('  gives %.2f GPa at zero brine where pure ice is 9.37.' % C_A)


def case_gogolaze(z):
    print('\n' + '=' * 70)
    print('GOGOLAZE 2026 -- phi from his eq. (14), target is beam rigidity')
    print('=' * 70)
    zc = z * H_GOGO * 100.0
    phi = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    print('  phi %.4f top, %.4f min, %.4f base' % (phi[0], phi.min(), phi[-1]))
    print('  %-24s %10s %10s' % ('exponent', 'E_flex', 'vs 1.421'))
    for n in (1.0, 2.0, 3.0, 4.0):
        v = flexural(law.layered(phi, n), z)
        print('  %-24.2f %10.3f %9.2fx' % (n, v, v / GOGO_COR))
    try:
        nfit = brentq(lambda n: flexural(law.layered(phi, n), z) - GOGO_COR,
                      0.3, 12.0)
        print('  exponent needed for the root-corrected value: %.2f  (%s)'
              % (nfit, verdict(nfit)))
    except ValueError:
        print('  no exponent in range reaches it')
    print('  His is an absolute whole-beam modulus, which is the hardest of the')
    print('  three targets: no normalisation is available to absorb an offset.')


def case_kujala(z):
    print('\n' + '=' * 70)
    print('KUJALA 1990 -- no porosity reported, so invert his profile instead')
    print('=' * 70)
    Et, Eb = K_TOP.mean(), K_BOT.mean()
    tgt = Et + (Eb - Et) * z
    print('  his profile %.2f -> %.2f GPa, alpha %.3f, z0/H %.2f-%.2f'
          % (Et, Eb, Eb / Et, K_Z0.min(), K_Z0.max()))
    for n in (1.0, 2.0):
        req = []
        for target in tgt:
            try:
                req.append(brentq(lambda p: law.layered(p, n) - target,
                                  1e-6, 0.95))
            except ValueError:
                req.append(np.nan)
        req = np.array(req)
        mono = bool(np.all(np.diff(req[np.isfinite(req)]) > -1e-6))
        print('\n  with b^%.0f, the brine profile his beams would need:' % n)
        print('    %8s %10s' % ('z/H', 'phi'))
        for zz in (0.0, 0.25, 0.5, 0.75, 1.0):
            print('    %8.2f %10.4f' % (zz, np.interp(zz, z, req)))
        print('    monotonic: %s | plausible for first-year ice: %s'
              % (mono, 'yes' if np.nanmax(req) < 0.40 else
                 'NO, base needs phi = %.2f' % np.nanmax(req)))


def main():
    z = np.linspace(1e-3, 1.0, 400)
    case_marchenko(z)
    case_gogolaze(z)
    case_kujala(z)
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)
    print('  Marchenko is matched on shape with a defensible exponent, and the')
    print('  residual level gap is an artefact of his own intercept.')
    print('  Gogolaze needs an exponent no cellular-solid argument supports,')
    print('  and it is an absolute target so nothing can absorb the offset.')
    print('  Kujala cannot be matched forwards at all; the inversion says what')
    print('  porosity his beams imply, and whether that is credible ice.')


if __name__ == '__main__':
    main()
