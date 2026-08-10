# -*- coding: utf-8 -*-
"""How strong are the single-packing claims really?

Several sweeps carry one cell per condition. Their effects can still be judged,
by transferring the packing scatter measured on the n=10 ensembles at comparable
microstructure, and by turning the n=4 nulls into bounds rather than assertions.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

os.chdir('C:/Users/stirpeg2/SpaX/SpaX/.claude/worktrees/skeletal-eringen-weibull/results')


def num(d, *c):
    for x in c:
        d[x] = pd.to_numeric(d[x], errors='coerce')
    return d


# ---- reference scatter on E_z/E_x, from the ten-packing column ensemble
ens = num(pd.read_csv('results_colseeds_all.csv'), 'E_x', 'E_z')
ens['slice'] = ens.run_id.str.extract(r'(z\d\d)')[0]
ens['a'] = ens.E_z / ens.E_x
sd = ens.groupby('slice')['a'].std(ddof=0)
free = sd.loc[['z05', 'z15', 'z25', 'z35', 'z45', 'z55', 'z65', 'z75']]
chan = sd.loc[['z85', 'z95']]
sd_free = float(free.mean())
sd_chan = float(chan.mean())
print('packing scatter on E_z/E_x, from n=10 ensembles')
print('   channel-free slices : %.4f  (range %.4f-%.4f)' % (sd_free, free.min(), free.max()))
print('   channelled slices   : %.4f  (range %.4f-%.4f)' % (sd_chan, chan.min(), chan.max()))
print()

# ---- single-packing sweeps, judged against that scatter
print('single-packing sweeps, effect size in units of the transferred scatter')
tests = []

m = num(pd.read_csv('results_morph.csv'), 'E_x', 'E_z')
m['a'] = m.E_z / m.E_x
tests.append(('morphology (sphericity)', m.a.max() - m.a.min(), sd_free, len(m)))

o = num(pd.read_csv('results_orient.csv'), 'E_x', 'E_z')
o['a'] = o.E_z / o.E_x
tests.append(('orientation (axis)', o.a.max() - o.a.min(), sd_free, len(o)))

p = num(pd.read_csv('results_perc.csv'), 'E_x', 'E_z', 'channel_vof_target')
p['a'] = p.E_z / p.E_x
poff = p[p.channel_vof_target == 0]
pon = p[p.channel_vof_target > 0]
tests.append(('percolation, channels off', poff.a.max() - poff.a.min(), sd_free, len(poff)))
tests.append(('percolation, channels on', pon.a.max() - pon.a.min(), sd_chan, len(pon)))

c = num(pd.read_csv('results_channel.csv'), 'E_x', 'E_z')
c['a'] = c.E_z / c.E_x
tests.append(('channel geometry', c.a.max() - c.a.min(), sd_chan, len(c)))

for name, eff, s, n in tests:
    print('   %-26s span %.4f = %5.1f sd   (n=%d cells, 1 each)'
          % (name, eff, eff / s, n))

# ---- the tilt null, expressed as a bound
print()
print('tilt sweep: what the null actually bounds')
t = {k: num(pd.read_csv('results_tilt%s.csv' % k), 'E_x', 'E_z') for k in ('00', '15', '30')}
r = {k: (v.E_z / v.E_x).dropna() for k, v in t.items()}
a, b = r['00'], r['30']
diff = a.mean() - b.mean()
se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
tc = stats.t.ppf(0.975, len(a) + len(b) - 2)
print('   straight %.4f (n=%d), 30 deg %.4f (n=%d)' % (a.mean(), len(a), b.mean(), len(b)))
print('   difference %+.4f, 95%% CI [%+.4f, %+.4f]' % (diff, diff - tc * se, diff + tc * se))
excess = a.mean() - 1.0
print('   the excess E_z/E_x - 1 is %.4f, so the CI bounds any dilution at %.0f%% of it'
      % (excess, 100 * (diff + tc * se) / excess))
