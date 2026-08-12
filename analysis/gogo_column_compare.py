import csv, re, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/analysis')
import ez_closure as ez

P = (r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull'
     r'/results/results_gogo_column.csv')
g = defaultdict(list)
for r in csv.DictReader(open(P, encoding='utf8', errors='replace')):
    m = re.match(r'GOGO_z(\d+)_s(\d+)', r['run_id'])
    if not m:
        continue
    try:
        ex, ph = float(r['E_x']), float(r['phi_inclusion'])
    except (ValueError, TypeError):
        continue
    if ex <= 0:
        continue
    g[int(m.group(1)) / 1000.0].append((ex / 1e9, ph))

z = np.array(sorted(g))
E_rve = np.array([np.mean([v[0] for v in g[k]]) for k in z])
phi = np.array([np.mean([v[1] for v in g[k]]) for k in z])
h = 1.0 / len(z)


def flex(E):
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    D = float(np.sum(E * h * ((z - z0) ** 2 + h ** 2 / 12.0)))
    return 12 * D, z0


print('all slices above phi_c = %.2f : %s' % (ez.PHI_C, bool((phi > ez.PHI_C).all())))
print('max phi = %.4f  (phi_0 = %.2f)' % (phi.max(), ez.PHI_0))
print()

E_ramp = np.array([float(ez.E_of_phi(p)) for p in phi])
Ep = ez.E_ICE * (1 - 1.65 * phi)
b = np.clip(1 - np.sqrt(np.clip(phi, 0, ez.PHI_0) / ez.PHI_0), 0, 1)
ne = ez.N_MID * (ez.A0_REF_MM / ez.A0_MM) ** ez.SPACING_EXP
E_step = np.maximum(np.where(phi >= ez.PHI_C, Ep * b ** ne, Ep), ez.E_FLOOR)

print('%-8s %8s %10s %10s %10s' % ('z/H', 'phi', 'E_rve', 'E_ramp', 'E_step'))
for i in range(len(z)):
    print('%-8.3f %8.4f %10.3f %10.3f %10.3f'
          % (z[i], phi[i], E_rve[i], E_ramp[i], E_step[i]))

print()
for lab, E in (('RVE cells (particle)', E_rve), ('closure, ramp', E_ramp),
               ('closure, step', E_step)):
    f, z0 = flex(E)
    print('%-22s E_flex=%7.3f GPa  z0/H=%.3f  alpha=%.3f'
          % (lab, f, z0, E[-1] / E[0]))

print('\nGogolaze measured : 0.785 - 1.421 GPa')
for lab, E in (('RVE cells', E_rve), ('ramp', E_ramp), ('step', E_step)):
    f, _ = flex(E)
    print('  %-10s %7.3f GPa -> %.1fx-%.1fx his   |  x0.49 -> %.3f (%.1fx-%.1fx)'
          % (lab, f, f / 1.421, f / 0.785, f * 0.49,
             f * 0.49 / 1.421, f * 0.49 / 0.785))
