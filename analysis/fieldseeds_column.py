"""fieldseeds and gogo_chanfrac, post-fix.

fieldseeds is the 15-packing-per-depth ensemble the field comparison of
section 4.6.2 rests on -- it is what makes the Marchenko level offset and the
matrix factor statements ensemble quantities rather than single-cell ones.

gogo_chanfrac sweeps how much of Gogolaze's brine is routed into channels.
"""
import csv, re, sys
from collections import defaultdict
import numpy as np

R = r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results/'


def load(name, pat):
    g = defaultdict(list)
    for r in csv.DictReader(open(R + name, encoding='utf8', errors='replace')):
        m = re.match(pat, r['run_id'])
        if not m:
            continue
        try:
            ex, ez, ph = float(r['E_x']), float(r['E_z']), float(r['phi_inclusion'])
        except (ValueError, TypeError):
            continue
        if ex <= 0:
            continue
        g[m.group(1)].append((ex / 1e9, ez / 1e9, ph))
    return g


print('=== fieldseeds: depth ensemble ===')
g = load('results_fieldseeds.csv', r'MSEED_z(\d+)_s\d+')
z = sorted(g, key=int)
print('%-6s %4s %18s %18s %9s' % ('z/H', 'n', 'E_x [GPa]', 'Ez/Ex', 'phi'))
E, S, Z, PH = [], [], [], []
for k in z:
    a = np.array(g[k])
    zc = int(k) / 100.0
    E.append(a[:, 0].mean()); S.append(a[:, 0].std())
    Z.append(zc); PH.append(a[:, 2].mean())
    print('%-6.2f %4d  %6.3f +- %.3f  %6.4f +- %.4f  %7.4f'
          % (zc, len(a), a[:, 0].mean(), a[:, 0].std(),
             (a[:, 1] / a[:, 0]).mean(), (a[:, 1] / a[:, 0]).std(),
             a[:, 2].mean()))
E = np.array(E); S = np.array(S); Z = np.array(Z)
print('\nCoV by depth (%%): ', np.round(100 * S / E, 2))
print('E_top %.3f  E_base %.3f  alpha %.4f' % (E[0], E[-1], E[-1] / E[0]))

h = 1.0 / len(E)
z0 = float(np.sum(E * h * Z) / np.sum(E * h))
D = float(np.sum(E * h * ((Z - z0) ** 2 + h ** 2 / 12.0)))
print('E_flex %.3f GPa   E_ext %.3f GPa   z0/H %.4f'
      % (12 * D, float(np.sum(E * h)), z0))
print('E_flex with x0.49 matrix: %.3f GPa   (Kujala measured 4.1 +- 0.5)'
      % (12 * D * 0.49))

print('\n=== gogo_chanfrac: brine routed into channels ===')
g2 = load('results_gogo_chanfrac.csv', r'GOCH_f(\d+)_s\d+')
print('%-8s %4s %18s %18s %9s' % ('frac', 'n', 'E_x [GPa]', 'Ez/Ex', 'phi'))
for k in sorted(g2, key=int):
    a = np.array(g2[k])
    print('%-8.2f %4d  %6.3f +- %.3f  %6.4f +- %.4f  %7.4f'
          % (int(k) / 100.0, len(a), a[:, 0].mean(), a[:, 0].std(),
             (a[:, 1] / a[:, 0]).mean(), (a[:, 1] / a[:, 0]).std(),
             a[:, 2].mean()))
