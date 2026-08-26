"""Which re-run results are compromised by having no slab primitive?

Every depth-column campaign meshes pockets and channels at every depth. Below
the percolation threshold that is the morphology section 4.4 shows to fail, so
any slice with phi > phi_c is over-predicted, and any column-level quantity
inherits that error weighted by how much the base contributes.

The question that matters most: the seasonal sweep. As a sheet warms, MORE
slices cross phi_c, so more of the column should switch to the layered
description. A particle-only column cannot show that, which means the
conclusion "warming does not deepen the critical layer" may be an artefact of
the missing primitive rather than a result.
"""
import csv, re, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/analysis')
import ez_closure as ez

R = r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results/'


def column(fname, pat):
    """-> z, E_particle, phi  (means over packings, sorted by depth)"""
    g = defaultdict(list)
    for r in csv.DictReader(open(R + fname, encoding='utf8', errors='replace')):
        m = re.match(pat, r['run_id'])
        if not m:
            continue
        try:
            ex, ph = float(r['E_x']), float(r['phi_inclusion'])
        except (ValueError, TypeError):
            continue
        if ex <= 0:
            continue
        g[int(m.group(len(m.groups())))].append((ex / 1e9, ph))
    ks = sorted(g)
    z = np.linspace(0, 1, len(ks) + 1)
    z = 0.5 * (z[:-1] + z[1:])
    E = np.array([np.mean([v[0] for v in g[k]]) for k in ks])
    P = np.array([np.mean([v[1] for v in g[k]]) for k in ks])
    return z, E, P


def hybrid(z, E, P):
    """Particle cells above phi_c, closure below -- the paper's own recipe."""
    Ec = np.array([float(ez.E_of_phi(p)) for p in P])
    return np.where(P >= ez.PHI_C, Ec, E)


def alpha_z0(z, E):
    h = 1.0 / len(E)
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    return E[-1] / E[0], z0


CASES = [
    ('colseeds  (FY, -20C)', 'results_colseeds.csv', r'CSEED_z(\d+)_s\d+'),
    ('seasonal  w20', 'results_seas.csv', r'SEAS_w20_z(\d+)'),
    ('seasonal  w12', 'results_seas.csv', r'SEAS_w12_z(\d+)'),
    ('seasonal  w06', 'results_seas.csv', r'SEAS_w06_z(\d+)'),
    ('fymy      FY', 'results_fymy.csv', r'FYMY_fy_z(\d+)'),
    ('fymy      MY', 'results_fymy.csv', r'FYMY_my_z(\d+)'),
    ('gogo_column', 'results_gogo_column.csv', r'GOGO_z(\d+)_s\d+'),
]

print('%-22s %5s %8s %8s %9s %9s %9s %9s' %
      ('campaign', 'n', 'phi_max', 'perc.', 'a_part', 'a_hyb', 'z0_part', 'z0_hyb'))
for lab, f, pat in CASES:
    try:
        z, E, P = column(f, pat)
    except Exception as e:
        print('%-22s  (%s)' % (lab, e))
        continue
    if len(E) == 0:
        print('%-22s  (no rows)' % lab)
        continue
    H = hybrid(z, E, P)
    ap, z0p = alpha_z0(z, E)
    ah, z0h = alpha_z0(z, H)
    nperc = int((P >= ez.PHI_C).sum())
    print('%-22s %5d %8.4f %6d/%-2d %9.4f %9.4f %9.4f %9.4f'
          % (lab, len(E), P.max(), nperc, len(P), ap, ah, z0p, z0h))

print('\nThe seasonal question: does the percolated fraction grow as it warms?')
for lab, f, pat in CASES[1:4]:
    z, E, P = column(f, pat)
    frac = (P >= ez.PHI_C).mean()
    H = hybrid(z, E, P)
    print('  %-14s percolated fraction of thickness = %4.0f%%   '
          'alpha particle %.3f -> hybrid %.3f'
          % (lab, 100 * frac, alpha_z0(z, E)[0], alpha_z0(z, H)[0]))
