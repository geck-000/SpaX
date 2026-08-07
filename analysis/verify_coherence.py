# -*- coding: utf-8 -*-
"""Whole-paper coherence audit for main_fix.tex.

Two independent checks:
  (a) VALUE  -- each canonical fact is recomputed from its source CSV and the
                paper's figure compared against it.
  (b) ECHO   -- each canonical fact is searched for across the whole .tex, and
                any *stale* variant of it that still appears is flagged. This is
                the check that catches a number updated in one section and left
                behind in the abstract, discussion or conclusions.
"""
import io, re, sys, os, glob
import numpy as np
import pandas as pd

TEX = 'C:/Users/stirpeg2/AppData/Local/Temp/overleaf-68d39c9d6e301aadbb376c0e/main_fix.tex'
RES = 'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results'
os.chdir(RES)
tex = io.open(TEX, encoding='utf8').read()


def num(df, col):
    return pd.to_numeric(df[col], errors='coerce')


# ---------------------------------------------------------------- recompute
ens = pd.read_csv('results_colseeds_all.csv')
ens['slice'] = ens['run_id'].str.extract(r'(z\d\d)')
ens['Ex'] = num(ens, 'E_x') / 1e9
ens['Ez'] = num(ens, 'E_z') / 1e9
g = ens.groupby('slice')
Ex = g['Ex'].mean()
ani = g.apply(lambda v: (v.Ez / v.Ex).mean(), include_groups=False)
anisd = g.apply(lambda v: (v.Ez / v.Ex).std(ddof=0), include_groups=False)
cov = g.apply(lambda v: 100 * v.Ex.std(ddof=0) / v.Ex.mean(), include_groups=False)

wb = pd.read_csv('results_weibull.csv').set_index('case')

perc = pd.read_csv('results_perc.csv')
tilt = {t: pd.read_csv('results_tilt%s.csv' % t) for t in ('00', '15', '30')}
tiltr = {t: (num(d, 'E_z') / num(d, 'E_x')).dropna() for t, d in tilt.items()}

bs = pd.concat([pd.read_csv(f) for f in
                ('results_basesweep_new.csv', 'results_basesweep_L065.csv',
                 'results_basesweep_L100.csv')], ignore_index=True)
for c in ('E_x', 'E_z', 'L'):
    bs[c] = num(bs, c)
bs = bs.dropna(subset=['E_x', 'E_z'])
bsg = bs.groupby('L')

er = pd.read_csv('results_eringen.csv')
erh = pd.read_csv('results_eringen_homog.csv')
eb = num(er, 'E_bending')
er = er[eb > 0]

ct = pd.read_csv('results_coltensor.csv')

# ---------------------------------------------------------------- facts
# (label, computed, printed-in-paper, list of stale strings that must be gone)
F = [
 ('surface modulus',      Ex['z05'],        8.92,  ['8.9~GPa', '$8.9$~GPa']),
 ('peak modulus',         Ex['z25'],        9.06,  []),
 ('base modulus',         Ex['z95'],        6.50,  ['4.85', '4.9$~GPa', r'\sim4.9']),
 ('aniso z05',            ani['z05'],       1.0043,['1.0005']),
 ('aniso z75',            ani['z75'],       1.0143,[]),
 ('aniso z85',            ani['z85'],       1.0505,['1.026\\pm0.004']),
 ('aniso base',           ani['z95'],       1.1324,['1.045\\pm0.011', '1.050\\pm0.012',
                                                    '1.059\\pm0.010']),
 ('aniso base sd',        anisd['z95'],     0.0129,[]),
 ('CoV surface',          cov['z05'],       0.21,  ['0.33\\%']),
 ('CoV base',             cov['z95'],       1.72,  ['1.50\\%']),
 ('tilt 0',               tiltr['00'].mean(),1.060,['1.041\\pm0.005']),
 ('tilt 15',              tiltr['15'].mean(),1.060,['1.029\\pm0.003']),
 ('tilt 30',              tiltr['30'].mean(),1.059,['1.031\\pm0.008']),
 ('base sweep L=0.50',    bsg['E_x'].mean()[0.50]/1e9, 6.58, ['4.76', '4.70']),
 ('base sweep aniso',     (bs[bs.L==0.50].E_z/bs[bs.L==0.50].E_x).mean(), 1.1198, []),
 ('SCF base P99',         wb.loc['BASE','P99_mean'], 3.59, ['5.60', '\\approx6.0',
                                                            'P99}\\approx5.6']),
 ('SCF ctrl P99',         wb.loc['CTRL','P99_mean'], 1.90, ['1.87', '1.82']),
 ('SCF base m=1',         wb.loc['BASE','SCFeff_m1'],1.25, ['1.66']),
 ('SCF base m=50',        wb.loc['BASE','SCFeff_m50'],4.11,['15.8']),
 ('E_bend widest',        er[er.L==0.80]['E_bending'].astype(float).mean()/1e9, 7.64, ['6.26']),
 ('E_bend homog widest',  num(erh,'E_bending')[num(erh,'L')==0.80].mean()/1e9, 11.30, []),
 ('coltensor base ratio', ct['E_ratio'].iloc[-1], 1.1135, []),
]

print('=' * 74)
print('(a) VALUE CHECK -- paper figure against recomputed source')
print('=' * 74)
print('%-22s %12s %12s   %s' % ('fact', 'computed', 'in paper', ''))
bad = 0
for lab, c, p, _ in F:
    c = float(c)
    tol = max(0.006 * abs(c), 0.006)
    ok = abs(c - p) <= tol
    if not ok:
        bad += 1
    print('%-22s %12.4f %12.4f   %s' % (lab, c, p, 'ok' if ok else '<-- MISMATCH'))

print()
print('=' * 74)
print('(b) ECHO CHECK -- stale variants still present anywhere in the .tex')
print('=' * 74)
echo = 0
for lab, c, p, stale in F:
    for s in stale:
        n = tex.count(s)
        if n:
            echo += 1
            print('  %-22s stale string %-22r still appears %dx' % (lab, s, n))
            for m in re.finditer(re.escape(s), tex):
                a = max(0, m.start() - 90)
                print('        ...%s...' % tex[a:m.end() + 60].replace('\n', ' '))
if not echo:
    print('  none -- no stale variant of any canonical fact remains')

print()
print('value mismatches: %d    stale echoes: %d' % (bad, echo))
