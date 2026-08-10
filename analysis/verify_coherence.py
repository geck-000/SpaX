# -*- coding: utf-8 -*-
"""Coherence audit, extended with the facts the first version did not cover:
comparisons (gas vs channels in the SCF table), the CLT neutral-plane numbers,
the cantilever fit, the Kujala beams, and the skeletal separation."""
import io, os, re
import numpy as np
import pandas as pd

TEX = 'C:/Users/stirpeg2/AppData/Local/Temp/overleaf-68d39c9d6e301aadbb376c0e/main_fix.tex'
RES = 'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results'
os.chdir(RES)
tex = io.open(TEX, encoding='utf8').read()


def num(d, *c):
    for x in c:
        d[x] = pd.to_numeric(d[x], errors='coerce')
    return d


ens = num(pd.read_csv('results_colseeds_all.csv'), 'E_x', 'E_z')
ens['slice'] = ens.run_id.str.extract(r'(z\d\d)')[0]
g = ens.groupby('slice')
Ex = g['E_x'].mean() / 1e9
ani = g.apply(lambda v: (v.E_z / v.E_x).mean(), include_groups=False)

wb = pd.read_csv('results_weibull.csv').set_index('case')
fail = num(pd.read_csv('results_failure.csv'), 'SCF_p99')
sk = num(pd.read_csv('results_skeletal.csv'), 'E_x', 'E_z', 'phi_soft_total')
sk = sk.dropna(subset=['E_x', 'E_z'])
sk['a'] = sk.E_z / sk.E_x
ct = pd.read_csv('results_coltensor.csv')

F = [
 ('surface modulus',        Ex['z05'],                    8.92),
 ('base modulus',           Ex['z95'],                    6.50),
 ('aniso surface',          ani['z05'],                   1.0043),
 ('aniso base',             ani['z95'],                   1.1324),
 ('SCF base P99 (weibull)', wb.loc['BASE', 'P99_mean'],   3.59),
 ('SCF ctrl P99',           wb.loc['CTRL', 'P99_mean'],   1.90),
 ('SCF gas P50',            wb.loc['GAS', 'SCFeff_m1'],   1.146),
 ('failure base P99',       fail.SCF_p99.iloc[-1],        3.57),
 ('failure top P99',        fail.SCF_p99.iloc[0],         2.024),
 ('skeletal aniso, cells',  sk.a.max(),                   1.79),
 ('skeletal phi max',       sk.phi_soft_total.max(),      0.45),
 ('coltensor base ratio',   ct.E_ratio.iloc[-1],          1.1135),
 ('coltensor base split',   100*(ct.E_y.iloc[-1]/ct.E_x.iloc[-1]-1), 0.74),
]

print('=' * 72)
print('(a) VALUE CHECK')
print('=' * 72)
bad = 0
for lab, c, p in F:
    c = float(c)
    ok = abs(c - p) <= max(0.006 * abs(c), 0.006)
    bad += (not ok)
    print('%-24s %11.4f %11.4f  %s' % (lab, c, p, 'ok' if ok else '<-- MISMATCH'))

print()
print('=' * 72)
print('(b) ECHO CHECK -- stale variants anywhere in the .tex')
print('=' * 72)
STALE = [
 ('4.85', 'old base modulus'), ('5.60', 'old SCF base'), ('20.0 &', 'old SCF max'),
 ('26.9', 'old SCF vol frac'), ('1.9\\,\\phi', 'old brine coeff'),
 ('2.2\\,\\phi', 'old gas coeff'), ('k=2.29', 'old cantilever k'),
 ('0.65\\,d', 'old couple-stress bound'), ('B/\\sqrt{AD}=0.12', 'old coupling'),
 ('3.3\\%$ predicted', 'old neutral offset'), ('falls at the top of the measured', 'kujala overstatement'),
 ('1.935', 'old SCF P90'), ('does not move at all', 'tilt overstatement'),
 ('$45\\%$ for the first-year', 'old FY drop'), ('$\\approx0.55$', 'old alpha'),
]
echo = 0
for pat, why in STALE:
    n = tex.count(pat)
    if n:
        echo += 1
        print('  %-26s %-34s %dx' % (why, repr(pat), n))
if not echo:
    print('  none')

print()
print('value mismatches: %d    stale echoes: %d' % (bad, echo))
