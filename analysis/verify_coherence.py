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



# The plate quantities and the Gogolaze forward evaluation. E_flex is the
# ABD-inverse modulus, 12(D - B^2/A), which is what Section 4.4's E_flex/E_ext
# and Section 4.5's E^t/E_flex both use; Table 6's "thickness average" is 12D
# about the geometric mid-plane instead, so both are checked rather than one
# being assumed to stand for the other.
def _column():
    import csv as _csv
    rows = list(_csv.DictReader(io.open('results_column_ensemble.csv', encoding='utf8')))
    Ex = np.array([float(r['E_x']) for r in rows]) / 1e9
    nu = np.array([float(r['nu_x']) for r in rows])
    h = 1.0 / len(rows)
    zc = (np.arange(len(rows)) + 0.5) * h
    T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
    S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])
    ph = _ez.brine_volume(T, S)
    Ec = np.array([float(_ez.E_of_phi(x)) for x in ph])
    E = np.where(ph > 0.05, Ec, Ex)
    Q = E / (1 - nu ** 2)
    A = np.sum(Q * h)
    B = np.sum(Q * h * (zc - 0.5))
    Dm = np.sum(Q * h * ((zc - 0.5) ** 2 + h ** 2 / 12.0))
    Mi = np.linalg.inv(np.array([[A, B], [B, Dm]]))
    return dict(Et=Ec[0], Eext=1 / Mi[0, 0], Eflex=12 / Mi[1, 1], D12mid=12 * Dm)


_plate = _column()


def _gogo(weight):
    """Beam rigidity of the Gogolaze column, 12D about its own neutral plane."""
    import csv as _csv
    import re as _re
    from collections import defaultdict as _dd
    g = _dd(list)
    for r in _csv.DictReader(io.open('results_gogo_column.csv', encoding='utf8',
                                     errors='replace')):
        m = _re.match(r'GOGO_z(\d+)_s(\d+)', r['run_id'])
        if not m:
            continue
        try:
            ex, ph = float(r['E_x']), float(r['phi_inclusion'])
        except (TypeError, ValueError):
            continue
        if ex > 0:
            g[int(m.group(1)) / 1000.0].append(ph)
    z = np.array(sorted(g))
    ph = np.array([np.mean(g[k]) for k in z])
    h = 1.0 / len(z)
    E = np.array([float(_ez.E_of_phi(x, weight=weight)) for x in ph])
    z0 = float(np.sum(E * h * z) / np.sum(E * h))
    return 12 * float(np.sum(E * h * ((z - z0) ** 2 + h ** 2 / 12.0)))


# The level correction, derived rather than asserted: c is the intersection of
# what the two whole-beam quantities demand, so a change in the closure or in
# either measured band moves it, and the manuscript's Eq. (7) has to move with
# it. That did not happen when the ramp was fixed, which is how 3.70 / 2.62
# survived in Table 6 after Table 7's ramp row had been rebuilt.
_c_kujala = (3.6 / _plate['D12mid'], 4.6 / _plate['D12mid'])
_c_gogo = (0.785 / _gogo('step'), 1.421 / _gogo('step'))

# Section 3.4's seasonal sweeps and Section 4.4's neutral-plane argument. The
# first two are closure-derived through the SEAS profiles, the third through the
# assembled column; all three had been left on a superseded configuration, the
# neutral-plane paragraph on the ramp while the paper adopts the step.
def _seas():
    import csv as _csv
    ph, Ec = {}, {}
    for r in _csv.DictReader(io.open('results_seas.csv', encoding='utf8')):
        for w in ('w20', 'w12', 'w06'):
            if 'SEAS_' + w in r['run_id']:
                ph.setdefault(w, []).append(float(r['phi_soft_total']))
                Ec.setdefault(w, []).append(float(r['E_x']) / 1e9)
    out = {}
    for w in ph:
        p_, c_ = np.array(ph[w]), np.array(Ec[w])
        E = np.array([float(_ez.E_of_phi(x)) for x in p_])
        Eh = np.where(p_ > 0.05, E, c_)
        out[w] = dict(alpha=Eh[-1] / c_[0], past=float((p_ > _ez.PHI_C).mean()))
    return out


_s = _seas()


def _z0(E, nu, h, zc):
    Q = E / (1 - nu ** 2)
    return float(np.sum(Q * h * zc) / np.sum(Q * h))


def _straight_line_z0():
    """Neutral plane of a straight line through the adopted column's endpoints."""
    import csv as _csv
    rows = list(_csv.DictReader(io.open('results_column_ensemble.csv', encoding='utf8')))
    Ex = np.array([float(r['E_x']) for r in rows]) / 1e9
    nu = np.array([float(r['nu_x']) for r in rows])
    h = 1.0 / len(rows)
    zc = (np.arange(len(rows)) + 0.5) * h
    T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
    S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])
    pp = _ez.brine_volume(T, S)
    Ec = np.array([float(_ez.E_of_phi(x)) for x in pp])
    E = np.where(pp > 0.05, Ec, Ex)
    return _z0(np.linspace(E[0], E[-1], len(E)), nu, h, zc)


# Section 3.5's mechanism numbers and the appendix sweeps, none of which had a
# check: the constriction and bracket exponents, the drainage ratios, the
# spacing exponent from the layer-count sweep, the above-phi_0 exponents, the
# anisotropy sweeps and the localisation percentiles.
def _pow(fn, key, state, xfun, col='E_x', scale=1e9):
    d = pd.read_csv(fn)
    d[col] = pd.to_numeric(d[col], errors='coerce') / scale
    d = d[d.run_id.str.contains(state)]
    x = d.run_id.map(xfun)
    m = d.groupby(x)[col].mean().dropna()
    lx, ly = np.log(m.index.values.astype(float)), np.log(m.values)
    return float(np.polyfit(lx, ly, 1)[0]), m


import re as _re
_nb = lambda r: float(_re.search(r'_n(\d+)', r).group(1))
_bb = lambda r: float(_re.search(r'_b(\d+)', r).group(1)) / 1000

_con, _ = _pow('results_bracket_nbridges.csv', 'n', '_drn', _nb)
_conu, _ = _pow('results_bracket_nbridges.csv', 'n', '_und', _nb)
_brk, _mbrk = _pow('results_bracket_bridge.csv', 'b', '_drn', _bb)
_brku, _mu = _pow('results_bracket_bridge.csv', 'b', '_und', _bb)
_nl = lambda r: 1.0 / float(_re.search(r'_n(\d+)', r).group(1))
_spc, _ = _pow('results_bracket_nlayers.csv', 'a0', '_drn', _nl)
_wl = pd.read_csv('results_weibull_layer_scf.csv').drop_duplicates('run_id')
_wl['g'] = _wl.run_id.str.replace(r'_s\d+$', '', regex=True)
_wg = _wl.groupby('g').SCF_p99.mean()
_wb = pd.read_csv('results_weibull.csv').set_index('case')

# Appendix A.1's replicate and size statistics, and the constituent sweeps. The
# thirteen layered conditions are the nine LSK cells and the four LCOL ones;
# nothing recorded that, and it is the only set of thirteen whose drained median
# and maximum are 1.10 and 3.87 and whose undrained pair is 5.55 and 19.2.
def _cov_set(files, state):
    out = []
    for fn in files:
        d = pd.read_csv(fn)
        d['E'] = pd.to_numeric(d['E_x'], errors='coerce')
        d['g'] = d.run_id.astype(str).str.replace(r'_s\d+$', '', regex=True)
        for g, x in d.groupby('g'):
            e = x.E.dropna()
            if len(e) >= 2 and state in g:
                out.append(100 * e.std(ddof=0) / e.mean())
    return np.array(sorted(out))


_LAY = ['results_layerskel.csv', 'results_layercol.csv']
_ld, _lu = _cov_set(_LAY, '_drn'), _cov_set(_LAY, '_und')


def _spread(fn, state, col='E_x'):
    d = pd.read_csv(fn)
    d[col] = pd.to_numeric(d[col], errors='coerce') / 1e9
    d['g'] = d.run_id.astype(str).str.replace(r'_s\d+$', '', regex=True)
    m = d[d.run_id.str.contains(state)].groupby('g')[col].mean()
    return 100 * m.std(ddof=0) / m.mean()


def _minspread(fn, tag):
    d = pd.read_csv(fn)
    d['E'] = pd.to_numeric(d['E_x'], errors='coerce') / 1e9
    q = d[d.run_id.str.contains(tag)]
    return 100 * (q.E.max() - q.E.min()) / q.E.min()


_bkc = pd.read_csv('results_brineKconst.csv')
_bkt = pd.read_csv('results_brineKtemp.csv')
for _x in (_bkc, _bkt):
    _x['E'] = pd.to_numeric(_x['E_x'], errors='coerce')

F += [
 ('layered CoV, drained median', float(np.median(_ld)),               1.10),
 ('layered CoV, drained max',   float(_ld.max()),                     3.87),
 ('layered CoV, undrained median', float(np.median(_lu)),             5.55),
 ('layered CoV, undrained max', float(_lu.max()),                     19.2),
 ('size, density held drained', _spread('results_bracket_density.csv', '_drn'), 1.868),
 ('size, density held undrained', _spread('results_bracket_density.csv', '_und'), 0.843),
 ('brine G sweep, pocket',      _minspread('results_brine.csv', 'iso_G'),   1.549),
 ('brine K sweep, pocket',      _minspread('results_brine.csv', 'iso_K'),   1.926),
 ('brine G sweep, channelled',  _minspread('results_brine.csv', 'chan_G'),  6.446),
 ('brine K sweep, channelled',  _minspread('results_brine.csv', 'chan_K'),  4.987),
 ('K(T) column shift %',        100 * (_bkt.E.mean() - _bkc.E.mean()) / _bkc.E.mean(), 0.0285),
 ('Goodier SCF at nu=0.33',     (27 - 15 * 0.33) / (2 * (7 - 5 * 0.33)),    2.06),
 ('constriction N^ (drained)', _con,                                 0.458),
 ('constriction N^ (undrained)', _conu,                              0.017),
 ('bracket b^ (drained)',    _brk,                                   0.689),
 ('bracket b^ (undrained)',  _brku,                                  0.007),
 ('drainage ratio at b=0.02', _mu.iloc[0] / _mbrk.iloc[0],           19.3),
 ('drainage ratio at b=0.28', _mu.iloc[-1] / _mbrk.iloc[-1],         2.9),
 ('spacing exponent a0^',    _spc,                                   0.69),
 ('SCF layered, phi=0.09',   _wg['WBLL_chan_drn'],                   6.23),
 ('SCF layered, phi=0.21',   _wg['WBLL_base_drn'],                   31.18),
 ('SCF layered sealed',      _wg['WBLL_base_und'],                   7.29),
 ('SCF lowest slice P99',    _wb.loc['BASE', 'P99_mean'],            3.593),
 ('SCF isolated pocket P99', _wb.loc['POCK', 'P99_mean'],            2.091),
 ('seas: alpha w20',        _s['w20']['alpha'],                      0.075),
 ('seas: alpha w12',        _s['w12']['alpha'],                      0.012),
 ('seas: alpha w06',        _s['w06']['alpha'],                      0.006),
 ('seas: past phi_c, w06',  _s['w06']['past'],                       0.30),
 ('neutral plane, straight', _straight_line_z0(),                    0.356),
 ('c: lower end (kujala)',  max(_c_kujala[0], _c_gogo[0]),           0.467),
 ('c: upper end (gogolaze)', min(_c_kujala[1], _c_gogo[1]),          0.552),
 ('c: gogolaze lower',      _c_gogo[0],                              0.305),
 ('c: kujala upper',        _c_kujala[1],                            0.596),
 ('sens: exponent +-2rms',  _E(n=_n - 2 * _ez.N_FIT_RMS) / _E(n=_n + 2 * _ez.N_FIT_RMS), 1.19),
 ('sens: a0 0.5-1.5 mm',    _E(a0_mm=1.5) / _E(a0_mm=0.5),           2.29),
 ('sens: phi_0 0.15-0.36',  _E(phi_0=0.36) / _E(phi_0=0.15),         5.38),
 ('sens: phi_0 0.14-0.23',  _E(phi_0=0.23) / _E(phi_0=0.14),         5.25),
 ('sens: b at phi=0.12',    _b,                                      0.225),
 ('N correction 4 -> 6',    _E(n_bridges=6) / _E(n_bridges=4),       1.22),
 ('N correction 4 -> 32',   _E(n_bridges=32) / _E(n_bridges=4),      2.81),
 ('plate: E_flex/E_ext',    _plate['Eflex'] / _plate['Eext'],        0.854),
 ('plate: E^t/E_flex',      _plate['Et'] / _plate['Eflex'],          1.214),
 ('plate: 12D about mid',   _plate['D12mid'],                        7.715),
 ('gogolaze rigidity, step', _gogo('step'),                          2.574),
 ('gogolaze rigidity, ramp', _gogo('ramp'),                          3.577),
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
 ('E_{\\mathrm{flex}}=1.22', 'pre-n(b) flexural ratio'),
 ('0.47$--$0.54', 'pre-rerun level correction'),
 ('0.47\\text{--}0.54', 'pre-rerun level correction in Eq. (7)'),
 ('3.70 / 2.62', 'pre-rerun gogolaze rigidity'),
 ('$0.14$--$0.24$', 'pre-rerun phi_0 bracket'),
 ('factor of two to three', 'rounded overshoot'),
 ('$\\times5.59$', 'superseded narrowed sensitivity'),
 ('$0.226$ to $0.089$', 'pre-rerun warming endpoint ratios'),
 ('$z_0/H=0.397$', 'neutral plane on the ramp column'),
 ('$0.06$ to $0.09$', 'ramp-column shortfall'),
 ('P99}=4.56', 'unlocated localisation percentile'),
 ('cell reaches $2.11$', 'peak quoted as a percentile'),
 ('$m\\approx4$', 'mis-stated m-norm crossover'),
 ('$1.030$--$1.037$', 'channels-off range on a channelled sweep'),
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
