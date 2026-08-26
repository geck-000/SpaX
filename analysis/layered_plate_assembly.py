"""What the plate looks like if the base is modelled with the morphology the
paper says belongs there.

The RVE column assembles ten pocket-and-channel cells. Section 4.4 argues the
particle description fails below the percolation depth, so the two warmest
slices are being represented by a morphology the paper itself rejects for them.
This substitutes the closure -- which carries the bridge factor below the
threshold and reduces to the pocket law above it -- and re-runs the neutral
plane.

Everything else is held: same ten slices, same thickness, same Poisson ratios,
same assembly.
"""
import csv, os, sys
import numpy as np

sys.path.insert(0, r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/analysis')
import ez_closure as ez

col = list(csv.DictReader(open(sys.argv[1], encoding='utf8')))
Ex_rve = np.array([float(r['E_x']) for r in col]) / 1e9
nu = np.array([float(r['nu_x']) for r in col])
n = len(col)
h = 1.0 / n
zc = (np.arange(n) + 0.5) * h

# the imposed (T,S) profile of Table 1
T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])
phi = ez.brine_volume(T, S)
b = 1.0 - np.sqrt(np.clip(phi, 0, ez.PHI_0) / ez.PHI_0)

E_clo = np.array([float(ez.E_of_phi(p)) for p in phi])

PHI_C = 0.05          # rule of fives: where the particle description ends


def z0_of(E):
    Q = E / (1.0 - nu ** 2)
    return float(np.sum(Q * h * zc) / np.sum(Q * h))


def assemble(E):
    Q = E / (1.0 - nu ** 2)
    A = np.sum(Q * h)
    z0 = z0_of(E)
    D = np.sum(Q * h * ((zc - z0) ** 2 + h ** 2 / 12.0))
    Bm = np.sum(Q * h * (zc - 0.5))
    Dm = np.sum(Q * h * ((zc - 0.5) ** 2 + h ** 2 / 12.0))
    return z0, A, 12 * D, abs(Bm) / np.sqrt(A * Dm), E[-1] / E[0]


print('%-5s %7s %6s %6s %9s %9s  %s' %
      ('z/H', 'phi_b', 'b', 'perc?', 'E_rve', 'E_closure', 'ratio'))
for i in range(n):
    print('%-5.2f %7.4f %6.3f %6s %9.3f %9.3f  %5.2f'
          % (zc[i], phi[i], b[i], 'YES' if phi[i] > PHI_C else '-',
             Ex_rve[i], E_clo[i], Ex_rve[i] / E_clo[i]))

# hybrid: keep the RVE cells where the particle description holds, use the
# closure only below the percolation depth, which is what the paper argues for
E_hyb = np.where(phi > PHI_C, E_clo, Ex_rve)

print('\n%-34s %8s %8s %8s %9s %8s' %
      ('column', 'z0/H', 'E_ext', 'E_flex', 'B/sqrtAD', 'alpha'))
for lab, E in (('RVE pockets+channels throughout', Ex_rve),
               ('closure throughout', E_clo),
               ('hybrid: closure below phi_c', E_hyb)):
    z0, A, Ef, cpl, al = assemble(E)
    print('%-34s %8.4f %8.3f %8.3f %9.4f %8.4f' % (lab, z0, A, Ef, cpl, al))

print('\nKujala measured: z0/H = 0.37-0.39, alpha = 0.12-0.19,')
print('                 E_top = 7.18-8.60 GPa, E_base = 0.86-1.56 GPa')
print('\nclosure base modulus  = %.3f GPa' % E_clo[-1])
print('closure surface       = %.3f GPa' % E_clo[0])
