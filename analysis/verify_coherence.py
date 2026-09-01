# -*- coding: utf-8 -*-
"""Coherence audit, extended with the facts the first version did not cover:
comparisons (gas vs channels in the SCF table), the CLT neutral-plane numbers,
the cantilever fit, the Kujala beams, and the skeletal separation."""
import io, os, re, sys
import numpy as np
import pandas as pd

# Paths were hardcoded to one workstation and to main_fix.tex, a manuscript that
# is no longer the one being submitted, so the audit silently checked the wrong
# file. Both are now arguments, defaulting to the current paper and this
# checkout's results:
#     SPAX_TEX=/path/to/new_rve_paper.tex python3 analysis/verify_coherence.py
TEX = os.environ.get(
    'SPAX_TEX',
    os.path.expanduser('~/opencode-workspace/sea-ice-paper/new_rve_paper.tex'))
RES = os.environ.get(
    'SPAX_RESULTS',
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
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
 ('coltensor base ratio',   ct.E_ratio.iloc[-1],          1.1135),
 ('coltensor base split',   100*(ct.E_y.iloc[-1]/ct.E_x.iloc[-1]-1), 0.74),
]
# The two skeletal checks that used to sit here asserted 1.79 and 0.45, numbers
# main_fix.tex carried and new_rve_paper.tex does not. A check against a claim
# the manuscript no longer makes reports a mismatch for ever, so they are gone
# rather than merely failing.

# --- the round-3 numbers ------------------------------------------------------
# Everything below is recomputed from the tensors and the closure rather than
# copied from the manuscript, so a change in either shows up as a mismatch.
import glob


def _C(fn):
    rows = []
    for line in io.open(fn, encoding='utf8'):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        q = line.split(',')
        if q[0] in ('11', '22', '33', '12', '13', '23'):
            rows.append([float(x) for x in q[1:7]])
        if len(rows) == 6:
            break
    return np.array(rows)


def _ratios(pattern):
    ry, rz = [], []
    for fn in sorted(glob.glob(pattern)):
        S = np.linalg.inv(_C(fn))
        Ex, Ey, Ez = 1 / S[0, 0], 1 / S[1, 1], 1 / S[2, 2]
        if Ex <= 0:
            continue
        ry.append(Ey / Ex)
        rz.append(Ez / Ex)
    return np.mean(ry), np.mean(rz)

T = os.path.join('..', 'tensors')
y85, z85 = _ratios(os.path.join(T, 'basetensor85_seeds', '*_z85_s*.csv'))
y95, z95 = _ratios(os.path.join(T, 'basetensor_seeds', '*_z95_s*.csv'))

# the layered basal cell, at fixed azimuth and averaged over azimuth
Cl = _C(os.path.join(T, 'layertensor', 'elasticity_tensor_LTEN_z95.csv'))
Sl = np.linalg.inv(Cl)
Sp = np.array([[Sl[0, 0], Sl[0, 1], 0], [Sl[0, 1], Sl[1, 1], 0], [0, 0, Sl[3, 3]]])
Ql = np.linalg.inv(Sp)
R = np.diag([1, 1, 2])
Qb = np.zeros((3, 3))
Sb = np.zeros((3, 3))
NTH = 7200
for th in np.linspace(0, 2 * np.pi, NTH, endpoint=False):
    c, sn = np.cos(th), np.sin(th)
    t = np.array([[c * c, sn * sn, 2 * c * sn],
                  [sn * sn, c * c, -2 * c * sn],
                  [-c * sn, c * sn, c * c - sn * sn]])
    ti = np.linalg.inv(t)
    Qb += ti @ Ql @ R @ t @ np.linalg.inv(R)
    Sb += R @ t @ np.linalg.inv(R) @ Sp @ ti
Qb /= NTH
Sb /= NTH
nur = -Sb[0, 1] / Sb[0, 0]
Q_reuss = (1 / Sb[0, 0]) / (1 - nur ** 2)

# The Section 5 sensitivity block. These were computed on the closure as it
# stood at 4281f3c, before n(b) replaced the single fitted exponent, and went
# unrevised for two rounds because nothing checked them: the sweep definitions
# live only in the prose, and 'a factor of three about the value the cells were
# solved at' turned out to mean 0.5-1.5 mm rather than 0.25-2.25. They are
# recomputed here so the next change to the closure cannot leave them behind.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import ez_closure as _ez                                             # noqa: E402

PHI_REF = 0.12
_b = 1 - np.sqrt(PHI_REF / _ez.PHI_0)
_n = _ez.n_of_b(_b)


def _E(**kw):
    return float(_ez.E_of_phi(PHI_REF, **kw))


F += [
 ('sens: exponent +-2rms',  _E(n=_n - 2 * _ez.N_FIT_RMS) / _E(n=_n + 2 * _ez.N_FIT_RMS), 1.19),
 ('sens: a0 0.5-1.5 mm',    _E(a0_mm=1.5) / _E(a0_mm=0.5),           2.29),
 ('sens: phi_0 0.15-0.36',  _E(phi_0=0.36) / _E(phi_0=0.15),         5.38),
 ('sens: phi_0 0.14-0.24',  _E(phi_0=0.24) / _E(phi_0=0.14),         5.59),
 ('sens: b at phi=0.12',    _b,                                      0.225),
 ('N correction 4 -> 6',    _E(n_bridges=6) / _E(n_bridges=4),       1.22),
 ('N correction 4 -> 32',   _E(n_bridges=32) / _E(n_bridges=4),      2.81),
 ('z85 in-plane E_y/E_x',   y85,                          0.998),
 ('z85 anisotropy E_z/E_x', z85,                          1.052),
 ('z95 in-plane E_y/E_x',   y95,                          0.999),
 ('z95 anisotropy E_z/E_x', z95,                          1.117),
 ('layered base E_y/E_x',   (1 / Sl[1, 1]) / (1 / Sl[0, 0]), 6.3),
 ('layered base nu_xy',     -Sl[0, 1] / Sl[0, 0],         0.047),
 ('layered base Q11 (GPa)', Ql[0, 0] / 1e9,               1.180),
 ('azimuthal Q, strain',    Qb[0, 0] / 1e9,               3.65),
 ('azimuthal Q, stress',    Q_reuss / 1e9,                1.84),
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
 ('B/\\sqrt{AD}=0.190', 'retired coupling'),
 ('upper half taking $66', 'old upper-half share'),
 ('0.999\\pm0.008', 'old z95 ensemble sd'),
 ('slightly\nmore steeply', 'reversed phase ordering'),
 ('\\times3.89', 'pre-n(b) phi_0 sensitivity'),
 ('\\times3.84', 'pre-n(b) narrowed phi_0 sensitivity'),
 ('\\times2.79', 'pre-n(b) spacing sensitivity'),
 ('\\times1.18$', 'pre-n(b) exponent sensitivity'),
 ('$0.93$--$1.04$', 'retired fixed exponent band'),
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
