"""Volume-weighted m-norm of the local stress field, layered against pocket.

A percentile is a numerical convenience; the statistic that carries a failure
probability is the Weibull integral, and its level depends strongly on m. The
paper already reports the sweep for the pocket-and-channel cells. This does the
same for the layered ones so the two are compared on the same footing rather
than on P99, which is where the descriptions were shown to separate.

    norm(m) = [ sum_e V_e (scf_e)^m / sum_e V_e ]^(1/m)

evaluated in log space: scf reaches ~114 here and 114^50 overflows a float.
Only tensile stress contributes to the Weibull integral, so the field is
clipped at zero before weighting.
"""
import glob
import os
import re
import collections

import numpy as np

DUMPS = 'weibull_dumps'
MS = (1, 2, 5, 10, 20, 50)


def mnorm(scf, vol, m):
    s = np.clip(scf, 0.0, None)
    good = s > 0
    if not good.any():
        return 0.0
    ls = np.log(s[good])
    lv = np.log(vol[good])
    a = m * ls + lv
    amax = a.max()
    lsum = amax + np.log(np.exp(a - amax).sum())
    return float(np.exp((lsum - np.log(vol.sum())) / m))


def main():
    groups = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(DUMPS, '*.npz'))):
        rid = os.path.basename(f)[:-4]
        cond = re.sub(r'_s\d+$', '', rid)
        groups[cond].append(f)

    print('%-20s %-4s %s' % ('condition', 'n',
                             '  '.join('m=%-6d' % m for m in MS)))
    for cond in sorted(groups):
        rows = []
        for f in groups[cond]:
            d = np.load(f)
            scf, vol = d['scf'], d['vol']
            rows.append([mnorm(scf, vol, m) for m in MS])
        a = np.array(rows)
        print('%-20s %-4d %s' % (cond, len(rows),
                                 '  '.join('%-8.2f' % v for v in a.mean(axis=0))))


if __name__ == '__main__':
    main()
