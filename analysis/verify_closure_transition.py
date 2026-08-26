"""Does the blended closure actually invert everywhere, and what does it do to
every number the paper quotes?"""
import csv, sys
import numpy as np
sys.path.insert(0, r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/analysis')
import ez_closure as ez

# --- 1. monotonicity, which is the reason for the change ------------------
phi = np.linspace(0.0, ez.PHI_0, 20001)
E = ez.E_of_phi(phi, floor=0.0)
d = np.diff(E)
print('MONOTONICITY over phi = 0 .. phi_0')
print('  strictly decreasing : %s' % bool(np.all(d < 0)))
if not np.all(d < 0):
    bad = np.where(d >= 0)[0]
    print('  non-decreasing at phi =', np.round(phi[bad][:6], 4))
print('  max jump between adjacent samples : %.5f GPa' % np.abs(d).max())
print('  E(0) = %.3f   E(phi_c) = %.3f   E(phi_0) = %.4f'
      % (E[0], float(ez.E_of_phi(ez.PHI_C, floor=0.0)), E[-1]))
print('  continuous at phi_c : %.6f vs %.6f'
      % (float(ez.E_of_phi(ez.PHI_C - 1e-9, floor=0.0)),
         float(ez.E_of_phi(ez.PHI_C + 1e-9, floor=0.0))))

# --- 2. the column and the plate ------------------------------------------
T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])
phic = ez.brine_volume(T, S)
col = list(csv.DictReader(open(r'C:/Users/stirpeg2/.claude/jobs/06dae8ab/tmp/column_new.csv')))
nu = np.array([float(r['nu_x']) for r in col])
Erve = np.array([float(r['E_x']) for r in col]) / 1e9
Ecl = np.array([float(ez.E_of_phi(p)) for p in phic])
Ehyb = np.where(phic >= ez.PHI_C, Ecl, Erve)
n = len(Ehyb); h = 1.0 / n; zc = (np.arange(n) + 0.5) * h


def assemble(Ev):
    Q = Ev / (1 - nu ** 2)
    A = np.sum(Q * h); z0 = float(np.sum(Q * h * zc) / A)
    D = np.sum(Q * h * ((zc - z0) ** 2 + h ** 2 / 12.0))
    Bm = np.sum(Q * h * (zc - 0.5))
    Dm = np.sum(Q * h * ((zc - 0.5) ** 2 + h ** 2 / 12.0))
    Ainv_ext = A
    return z0, Ainv_ext, 12 * D, abs(Bm) / np.sqrt(A * Dm), Ev[-1] / Ev[0]


print('\nCOLUMN (blend below phi_c, RVE cells above)')
print('  E per slice:', np.round(Ehyb, 3))
for lab, Ev in (('particle throughout', Erve), ('blend below phi_c', Ehyb)):
    z0, Aext, Ef, cpl, al = assemble(Ev)
    print('  %-22s E_base %6.3f  alpha %.4f  z0/H %.4f  E_flex %.3f  B/sqrtAD %.4f'
          % (lab, Ev[-1], al, z0, Ef, cpl))
lin = np.linspace(Ehyb[0], Ehyb[-1], n)
print('  straight line through the blend endpoints: z0/H = %.4f'
      % assemble(lin)[0])
print('  Kujala: E_base 0.86-1.56, alpha 0.12-0.19, z0/H 0.37-0.39')

# --- 3. the three field comparisons ---------------------------------------
z = np.linspace(1e-3, 1.0, 400)
tgt = 4.4 * (1 - 0.62 * z ** 0.6)
phi_m = (np.log(7.23 / tgt) / 4.2) ** 2
from scipy.optimize import minimize_scalar
best = minimize_scalar(lambda q: float(np.mean(
    (ez.E_of_phi(phi_m, n=q, floor=0.0) / float(ez.E_of_phi(phi_m[0], n=q, floor=0.0))
     - tgt / tgt[0]) ** 2)), bounds=(0.05, 8.0), method='bounded').x
print('\nMARCHENKO  his alpha 0.384')
for lab, q in (('n_lo', ez.N_LO), ('n_mid', ez.N_MID), ('n_hi', ez.N_HI),
               ('best %.2f' % best, best)):
    Em = ez.E_of_phi(phi_m, n=q, floor=0.0)
    print('  %-12s alpha %.3f   E_top %.2f' % (lab, Em[-1] / Em[0], Em[0]))

zc2 = z * 32.0
phi_g = (0.29315 * zc2 ** 2 - 5.124 * zc2 + 85.977) / 1000.0
print('\nGOGOLAZE  measured 0.785-1.421 GPa')
for lab, q in (('n_lo', ez.N_LO), ('n_mid', ez.N_MID), ('n_hi', ez.N_HI)):
    print('  %-6s beam rigidity %.3f GPa' % (lab, ez.flexural(
        ez.E_of_phi(phi_g, n=q, floor=0.0), z)))

print('\nKUJALA inversion: does every measured modulus now have a preimage?')
Emax = float(ez.E_of_phi(0.0, floor=0.0))
for t in (8.60, 8.05, 7.18, 4.0, 1.56, 1.16, 0.86):
    lo, hi = 1e-6, ez.PHI_0 - 1e-9
    f = lambda p: float(ez.E_of_phi(p, floor=0.0)) - t
    ok = f(lo) * f(hi) < 0
    print('  E = %5.2f GPa -> %s' % (t, 'solvable' if ok else 'NO SOLUTION'))
