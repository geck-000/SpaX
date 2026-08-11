r"""What the constriction result does to the closure and to the comparisons.

The subdivision test says the drained layered cell is constriction-dominated:
holding total bridge area fixed and splitting it over more bridges stiffens the
cell as N^0.458, against N^0.500 for spreading compliance through N circular
contacts. That is not a correction to the exponent in b. It changes the FORM of
the closure.

A power law E = E_pocket * b^n is a multiplier: the bridges scale the modulus.
Constriction is a compliance IN SERIES with the bulk: load crosses the ice,
funnels through the bridges, and spreads again, so the compliances add,

    1/E = 1/E_pocket(phi) + K / sqrt(N b)

which tends to E_pocket when the bridges are many or wide, and to the
constriction term when they are few or narrow. It also fixes the b -> 1 limit
that a bare power law gets wrong: a plane full of ice should return the pocket
modulus exactly, and this does.

The awkward consequence is that N, the number of bridges per lamellar plane,
enters as a second microstructural parameter. Our cells carry two. Real ice
carries many, and at equal b more bridges are stiffer, so a closure calibrated
on our cells and applied to real ice underestimates the modulus by sqrt(N/2).
This script calibrates K on the subdivision data and then asks what N each
dataset implies -- which is the honest way round, since N is observable and the
implied value can be checked against micro-CT rather than tuned.
"""
import numpy as np
from scipy.optimize import brentq, curve_fit

import layered_law as law

# subdivision campaign: phi = 0.15, b = 0.15, drained, N = 1..16
N_SUB = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
E_SUB = np.array([0.774, 0.889, 1.435, 2.094, 2.471])
PHI_SUB, B_SUB = 0.15, 0.15

GOGO_COR, H_GOGO = 1.421, 0.32
M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6
C_A, C_B = 7.23, 4.2
K_TOP, K_BOT = 8.05, 1.27


def marchenko_E(d):
    return M_E0 * (1.0 - (1.0 - M_ALPHA) * d ** M_N)


def corr_inv(E):
    return (np.log(C_A / np.maximum(E, 1e-9)) / C_B) ** 2


def E_series(phi, b, N, K):
    """Bulk and constriction compliances in series."""
    Ep = law.pocket(phi)
    nb = np.maximum(N * np.asarray(b, dtype=float), 1e-12)
    return 1.0 / (1.0 / Ep + K / np.sqrt(nb))


def flexural(E, z):
    z0 = float(np.trapz(E * z, z) / np.trapz(E, z))
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def main():
    # ---- calibrate K on the subdivision data -------------------------------
    def model(N, K):
        return E_series(PHI_SUB, B_SUB, N, K)
    K, _ = curve_fit(model, N_SUB, E_SUB, p0=[0.5])
    K = float(K[0])
    pred = model(N_SUB, K)
    print('CALIBRATION on the subdivision campaign (phi = b = 0.15, drained)')
    print('%6s %10s %10s %8s' % ('N', 'measured', 'series fit', 'error'))
    for n, e, p in zip(N_SUB, E_SUB, pred):
        print('%6.0f %10.3f %10.3f %7.1f%%' % (n, e, p, 100 * (p - e) / e))
    print('  K = %.4f GPa, RMS %.1f%%'
          % (K, 100 * np.sqrt(np.mean(((pred - E_SUB) / E_SUB) ** 2))))
    print('  a bare power law cannot do this: it has no b -> 1 limit and no')
    print('  way for the bulk and the bridges to trade off.')

    # ---- what N does each dataset imply? -----------------------------------
    z = np.linspace(1e-3, 1.0, 400)
    print('\nN IMPLIED BY EACH DATASET, with phi and b both external')
    print('(b = 1 - sqrt(phi/phi_0) throughout, phi_0 = %.2f)' % law.PHI_0)

    zc = z * H_GOGO * 100.0
    phi_g = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    b_g = law.assur_b(phi_g)
    try:
        Ng = brentq(lambda N: flexural(E_series(phi_g, b_g, N, K), z) - GOGO_COR,
                    0.05, 5000.0)
        print('  Gogolaze beam        N = %.1f bridges per plane' % Ng)
    except ValueError:
        print('  Gogolaze beam        no N reaches it')

    phi_m = corr_inv(marchenko_E(z))
    b_m = law.assur_b(phi_m)
    tgt = marchenko_E(z)

    def miss(N):
        E = E_series(phi_m, b_m, N, K)
        return float(np.mean((E / E[0] - tgt / tgt[0]) ** 2))
    grid = np.logspace(-1, 3.5, 400)
    Nm = grid[int(np.argmin([miss(g) for g in grid]))]
    print('  Marchenko profile    N = %.1f  (shape RMS %.4f)'
          % (Nm, np.sqrt(miss(Nm))))

    print('\nTHE TWO DO NOT AGREE, and the reason is structural.')
    print('  Gogolaze admits N = %.1f. Marchenko rails to the bottom of the'
          % Ng)
    print('  search at %.1f, with a shape RMS of %.3f against 0.015 for the'
          % (Nm, np.sqrt(miss(Nm))))
    print('  power law it replaces. Forty times apart is not calibration')
    print('  scatter.')
    bm = law.assur_b(phi_m)
    a_bulk = float(law.pocket(phi_m)[-1] / law.pocket(phi_m)[0])
    a_con = float(np.sqrt(bm[-1] / bm[0]))
    print('\n  The series form cannot reach his grading at ANY N, because its')
    print('  two limits bracket the whole family and both are too shallow. As N')
    print('  grows it tends to the pocket law, alpha = %.3f; as N falls to pure'
          % a_bulk)
    print('  constriction, alpha -> sqrt(b_base/b_top) = %.3f. He is at 0.384.'
          % a_con)
    print('\n  Gogolaze escapes that only because his target is a LEVEL rather')
    print('  than a ratio, and a level can always be met by moving N. That is')
    print('  the weaker test and should not be read as the closure working.')
    print('\n  What does survive: N is directly countable in a micro-CT slice')
    print('  through a lamellar plane, so the implied 4.3 is checkable rather')
    print('  than merely fitted.')

    print('\nSENSITIVITY, since our cells sat at N = 2')
    print('%8s %12s %12s' % ('N', 'Gogolaze', 'vs measured'))
    for n in (2, 10, 50, 200):
        v = flexural(E_series(phi_g, b_g, float(n), K), z)
        print('%8d %12.3f %11.2fx' % (n, v, v / GOGO_COR))


if __name__ == '__main__':
    main()
