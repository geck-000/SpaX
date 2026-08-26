r"""Rederive Table tab:nbridge: the bridge count against the field comparisons.

The three columns are the assembled reference column (Table 1 of the paper) run
through the closure at each count N, plus the Gogolaze beam-3 flexural modulus
after the level correction c:

    alpha   = (E_base/(1-nu^2)) / (E_top/(1-nu^2))   basal-to-surface ratio
    z0/H    = neutral plane of the Q = E/(1-nu^2) assembly
    cE_app  = c * E_flex(Gogolaze beam 3),  c = 0.54  (midpoint of 0.49-0.59)

The count enters through the (N/N_cells)^(kappa w) factor of Eq. (nbridge),
which the closure retains only to state what is not known. This prints the
table with the CURRENT closure so the paper's tab:nbridge can be re-derived.

    python3 analysis/nbridge_table.py
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ez_closure as ez

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

T_REF = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S_REF = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])

H_GOGO = 0.32                       # m, beam 3
C_LEVEL = 0.54                      # midpoint of the 0.49-0.59 band

COUNTS = (2, 3, 4, 5, 10, 30)


def ref_nu():
    """Poisson ratios of the reference-column RVE cells, in slice order."""
    path = os.path.join(ROOT, 'results', 'results_column.csv')
    return np.array([float(r['nu_x']) for r in csv.DictReader(open(path))])


def alpha_z0(E, nu):
    """Basal-to-surface ratio and neutral plane of the Q = E/(1-nu^2) assembly."""
    n = len(E)
    h = 1.0 / n
    z = (np.arange(n) + 0.5) * h
    Q = E / (1.0 - nu ** 2)
    z0 = float(np.sum(Q * h * z) / np.sum(Q * h))
    alpha = float(Q[-1] / Q[0])
    return alpha, z0


def gogo_flex(n_bridges):
    """Flexural modulus (GPa) the closure gives on the Gogolaze beam-3 profile."""
    z = np.linspace(0.001, 0.999, 600)
    zc = z * H_GOGO * 100.0
    phi = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    E = np.array([float(ez.E_of_phi(p, n_bridges=n_bridges)) for p in phi])
    z0 = np.trapz(E * z, z) / np.trapz(E, z)
    D = np.trapz(E * (z - z0) ** 2, z)
    return 12.0 * D / (z[-1] - z[0]) ** 3


def main():
    phi = ez.brine_volume(T_REF, S_REF)
    nu = ref_nu()

    print('reference column phi:', np.round(phi, 4))
    print()
    print('%-9s %10s %9s %12s' % ('N', 'alpha', 'z0/H', 'cE_app (GPa)'))
    print('-' * 44)
    for N in COUNTS:
        E = np.array([float(ez.E_of_phi(p, n_bridges=N)) for p in phi])
        alpha, z0 = alpha_z0(E, nu)
        cE = C_LEVEL * gogo_flex(N)
        tag = ' (cells)' if N == ez.N_CELLS else ''
        print('%-9s %10.3f %9.3f %12.2f%s'
              % (str(N) + tag, alpha, z0, cE, ''))


if __name__ == '__main__':
    main()
