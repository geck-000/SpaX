"""Re-read every fresh campaign the paper still quotes from pre-fix runs."""
import csv, os, re
from collections import defaultdict
import numpy as np

R = r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results/'


def rows(name):
    out = []
    p = R + name
    if not os.path.exists(p):
        return out
    for r in csv.DictReader(open(p, encoding='utf8', errors='replace')):
        try:
            ex, ez = float(r['E_x']), float(r['E_z'])
        except (ValueError, TypeError, KeyError):
            continue
        if ex <= 0:
            continue
        try:
            ph = float(r.get('phi_inclusion') or 'nan')
        except ValueError:
            ph = float('nan')
        out.append((r['run_id'], ex / 1e9, ez / 1e9, ph))
    return out


def grp(name, pat, keyf=lambda m: m.group(1)):
    g = defaultdict(list)
    for rid, ex, ez, ph in rows(name):
        m = re.match(pat, rid)
        if m:
            g[keyf(m)].append((ex, ez, ph))
    return g


def show(title, g, fmt='%s'):
    print('\n=== %s ===' % title)
    if not g:
        print('  (no rows matched)')
        return
    for k in sorted(g):
        a = np.array(g[k])
        r = a[:, 1] / a[:, 0]
        print('  %-14s n=%2d  E_x=%7.3f +- %.3f   Ez/Ex=%.4f +- %.4f  phi=%.4f'
              % (fmt % k, len(a), a[:, 0].mean(), a[:, 0].std(),
                 r.mean(), r.std(), np.nanmean(a[:, 2])))


# --- 4.2.1 seasonal + multi-year: the two remaining [RERUN] markers --------
show('SEASONAL (4.2.1 [RERUN])', grp('results_seas.csv', r'(SEAS_[A-Za-z0-9]+?)_z\d+'))
print('  --- per surface temperature, by depth ---')
g = grp('results_seas.csv', r'SEAS_(\w+?)_z(\d+)', lambda m: (m.group(1), int(m.group(2))))
tops = sorted({k[0] for k in g})
for t in tops:
    ks = sorted([k for k in g if k[0] == t], key=lambda k: k[1])
    if not ks:
        continue
    E = [np.mean([v[0] for v in g[k]]) for k in ks]
    print('    %-8s top %.3f  base %.3f  alpha %.3f  (%d depths)'
          % (t, E[0], E[-1], E[-1] / E[0], len(ks)))

show('FY vs MULTI-YEAR (4.2.1 [RERUN])', grp('results_fymy.csv', r'(FYMY_\w+?)_z\d+'))
g = grp('results_fymy.csv', r'FYMY_(\w+?)_z(\d+)', lambda m: (m.group(1), int(m.group(2))))
for t in sorted({k[0] for k in g}):
    ks = sorted([k for k in g if k[0] == t], key=lambda k: k[1])
    E = [np.mean([v[0] for v in g[k]]) for k in ks]
    if E:
        print('    %-8s top %.3f  base %.3f  drop %.1f%%'
              % (t, E[0], E[-1], 100 * (1 - E[-1] / E[0])))

# --- 4.2.2 tilt dilution ---------------------------------------------------
print('\n=== TILT DILUTION (4.2.2) ===')
for f, lab in (('results_tilt00.csv', '0 deg'), ('results_tilt15.csv', '15 deg'),
               ('results_tilt30.csv', '30 deg')):
    a = np.array([(x[1], x[2]) for x in rows(f)])
    if len(a):
        r = a[:, 1] / a[:, 0]
        print('  %-8s n=%d  E_x=%.3f +- %.3f   Ez/Ex=%.4f +- %.4f'
              % (lab, len(a), a[:, 0].mean(), a[:, 0].std(), r.mean(), r.std()))

# --- 4.2.3 constituents ----------------------------------------------------
show('GAS SWEEP (4.2.3)', grp('results_gas_dilute.csv', r'(GASD?_\w+)_s\d+'))
show('BRINE MODULI (4.2.3)', grp('results_brine.csv', r'(BR\w+?)_s?\d*$'))
for f, lab in (('results_brineKconst.csv', 'K const'),
               ('results_brineKtemp.csv', 'K(T)')):
    a = np.array([(x[1], x[2]) for x in rows(f)])
    if len(a):
        print('  %-10s n=%2d  mean E_x %.4f GPa' % (lab, len(a), a[:, 0].mean()))

# --- 4.1 cell size ---------------------------------------------------------
show('BASE CELL SIZE (4.1)', grp('results_basesweep_new.csv', r'BSW_L(\d+)'))
show('CHANNEL CELL SIZE (4.1)', grp('results_sizechan.csv', r'SZCH_L(\d+)'))

# --- 4.6 column variants ---------------------------------------------------
show('SKELETAL (4.6)', grp('results_skeletal.csv', r'(SKEL_\w+?)_?s?\d*$'))
show('STEEP COLUMN (4.6)', grp('results_steep_column.csv', r'STEEP_z(\d+)'))
show('LOWBASE (4.6)', grp('results_lowbase.csv', r'(LOWB_\w+?)_z\d+'))
show('SALINITY FAMILY (4.6.2)', grp('results_salfamily.csv', r'(SAL_\w+?)_z\d+'))
