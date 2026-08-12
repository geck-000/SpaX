"""What does Gogolaze's beam say about phi_0?

His base slice carries phi = 0.2084, just past Assur's phi_0 = 0.20, so b = 0
there and the closure returns its floor. But his beam demonstrably has finite
rigidity, so the question is whether phi_0 = 0.20 is consistent with that, and
what range of phi_0 is.

phi_0 is the paper's ASSUMED ingredient and its dominant uncertainty (x3.5 over
0.15-0.36). If his beam bounds it from below, that range narrows.
"""
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
phi = np.array([np.mean([v[1] for v in g[k]]) for k in z])
h = 1.0 / len(z)
BEAM = 0.49
MEAS = (0.785, 1.421)


def E_of(phi, phi_0, mode):
    Ep = ez.E_ICE * (1 - 1.65 * phi)
    b = np.clip(1 - np.sqrt(np.clip(phi, 0, phi_0) / phi_0), 0, 1)
    ne = ez.N_MID * (ez.A0_REF_MM / ez.A0_MM) ** ez.SPACING_EXP
    if mode == 'ramp':
        w = np.clip((phi - ez.PHI_C) / (phi_0 - ez.PHI_C), 0, 1)
        E = Ep * b ** (ne * w)
    else:
        E = np.where(phi >= ez.PHI_C, Ep * b ** ne, Ep)
    return np.maximum(E, ez.E_FLOOR)


def flex(E):
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    return 12 * float(np.sum(E * h * ((z - z0) ** 2 + h ** 2 / 12.0)))


print('his base slice phi = %.4f    Assur phi_0 = %.2f' % (phi.max(), ez.PHI_0))
print('slices at or above phi_0=0.20 : %d of %d'
      % (int((phi >= 0.20).sum()), len(phi)))
print()
print('%-8s %10s %10s %12s %12s'
      % ('phi_0', 'b(base)', 'E_base', 'Eflex ramp', 'Eflex step'))
for p0 in (0.15, 0.20, 0.21, 0.25, 0.30, 0.36):
    bb = max(0.0, 1 - np.sqrt(min(phi.max(), p0) / p0))
    Er = E_of(phi, p0, 'ramp'); Es = E_of(phi, p0, 'step')
    print('%-8.2f %10.4f %10.3f %12.3f %12.3f'
          % (p0, bb, Er[-1], flex(Er) * BEAM, flex(Es) * BEAM))

print('\nhis measured beam rigidity: %.3f - %.3f GPa' % MEAS)
print('(all values above already carry the %.2f beam-effective matrix factor)'
      % BEAM)

# which phi_0 put the ramp inside his band?
print('\nphi_0 consistent with his band, ramp form:')
ok = []
for p0 in np.arange(0.16, 0.40, 0.005):
    f = flex(E_of(phi, p0, 'ramp')) * BEAM
    if MEAS[0] <= f <= MEAS[1]:
        ok.append(p0)
print('   %.3f - %.3f' % (min(ok), max(ok)) if ok else '   none in 0.16-0.40')
ok2 = []
for p0 in np.arange(0.16, 0.40, 0.005):
    f = flex(E_of(phi, p0, 'step')) * BEAM
    if MEAS[0] <= f <= MEAS[1]:
        ok2.append(p0)
print('phi_0 consistent with his band, step form:')
print('   %.3f - %.3f' % (min(ok2), max(ok2)) if ok2 else '   none in 0.16-0.40')

# how much does the base slice actually matter to the rigidity?
print('\nhow much does the base slice contribute to D?')
for mode in ('ramp', 'step'):
    E = E_of(phi, ez.PHI_0, mode)
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    contrib = E * h * ((z - z0) ** 2 + h ** 2 / 12.0)
    print('   %-5s base slice = %.2f%% of D;  drop it entirely -> Eflex %.3f'
          % (mode, 100 * contrib[-1] / contrib.sum(),
             12 * contrib[:-1].sum() * BEAM))
