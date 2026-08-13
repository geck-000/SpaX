"""The bending size-effect measurement, re-read post-fix.

Section 4.3.1 regresses the apparent bending modulus on 1/L^2, the form both
nonclassical families take, and refers it to a matched set of homogeneous cells
which carry no microstructure and so measure the bias of the extraction itself.
"""
import csv, re
from collections import defaultdict
import numpy as np

R = r'C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results/'


def load(f, pat):
    g = defaultdict(list)
    for r in csv.DictReader(open(R + f, encoding='utf8', errors='replace')):
        m = re.match(pat, r['run_id'])
        if not m:
            continue
        try:
            e = float(r['E_bending'])
        except (ValueError, TypeError, KeyError):
            continue
        if e > 0:
            g[int(m.group(1))].append(e / 1e9)
    return g


ch = load('results_eringen.csv', r'ERG_L(\d+)_s\d+')
ho = load('results_eringen_homog.csv', r'ERGH_L(\d+)')
Ls = sorted(set(ch) & set(ho))
d = 0.06                      # mean inclusion diameter, model units

print('%-8s %-8s %5s %16s %12s %12s'
      % ('L', 'L/d', 'n', 'channelled', 'homogeneous', 'corrected'))
Lm, Ec, Eh, Ecorr, allpts = [], [], [], [], []
for L in Ls:
    a = np.array(ch[L]); h = np.mean(ho[L])
    Lv = L / 1000.0
    corr = a.mean() / h
    Lm.append(Lv); Ec.append(a.mean()); Eh.append(h); Ecorr.append(corr)
    for v in a:
        allpts.append((1.0 / Lv ** 2, v / h))
    print('%-8.3f %-8.1f %5d  %6.3f +- %.3f  %10.3f  %10.4f'
          % (Lv, Lv / d, len(a), a.mean(), a.std(), h, corr))

Lm = np.array(Lm); Ec = np.array(Ec); Eh = np.array(Eh)
print('\nhomogeneous control drifts %.1f%% across the sweep'
      % (100 * (Eh.max() - Eh.min()) / Eh.min()))
print('channelled cells drift     %.1f%%'
      % (100 * (Ec.max() - Ec.min()) / Ec.min()))

x = np.array([p[0] for p in allpts]); y = np.array([p[1] for p in allpts])
n = len(x)
p, cov = np.polyfit(x, y, 1, cov=True)
se = np.sqrt(cov[0, 0])
tstat = p[0] / se
from math import erf, sqrt
pval = 2 * (1 - 0.5 * (1 + erf(abs(tstat) / sqrt(2))))
print('\ncorrected regression on 1/L^2, all %d packings' % n)
print('  slope = %+.4f +- %.4f   (t = %.2f, p = %.3f)' % (p[0], se, tstat, pval))
print('  95%% CI = [%+.4f, %+.4f]' % (p[0] - 1.96 * se, p[0] + 1.96 * se))

# bounds. couple-stress: slope = 12 mu l^2 ; Eringen: slope = -E_inf (e0a)^2
mu, Einf = 1.97, np.mean(Ec)
hi = p[0] + 1.96 * se
lo = p[0] - 1.96 * se
if hi > 0:
    print('  couple-stress   l    < %.3f d' % (np.sqrt(hi * Einf / (12 * mu)) / d))
if lo < 0:
    print('  Eringen         e0a  < %.3f d' % (np.sqrt(-lo) / d))
