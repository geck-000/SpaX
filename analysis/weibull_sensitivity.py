"""Weakest-link sensitivity of the localisation measure (plain python3, no Abaqus).

Reads the per-element (scf, volume) dumps written by scf_extract.py and reports
how the effective stress-concentration factor depends on the Weibull modulus m,
so that the P99 quoted in the manuscript can be placed on a physical scale
rather than defended as a numerical convenience.

The measure
-----------
For a weakest-link solid the failure probability of a body under a field
sigma_1(x) is

    P_f = 1 - exp[ -(1/V0) INT (sigma_1/sigma_0)^m dV ]

so the field enters only through the m-norm of the normalised stress. Defining
an effective concentration factor as the volume-weighted m-norm of the SCF,

    SCF_eff(m) = [ (1/V) INT SCF^m dV ]^(1/m)

gives a single number that is comparable across cases, reduces to the
volume-weighted mean at m=1, and tends to the absolute maximum as m -> inf.
Each percentile of the SCF distribution equals SCF_eff at some particular m, so
reporting SCF_eff(m) makes explicit the choice that quoting P99 makes silently.

Usage
-----
    python3 weibull_sensitivity.py <dump_dir> [out_prefix]

<dump_dir> holds the .npz files from `scf_extract.py ... dump.npz`, named
WBL_<CASE>_s<seed>.npz. Writes <out_prefix>.csv and <out_prefix>.png
(default: results_weibull).
"""
import glob
import os
import re
import sys

import numpy as np

# m from mild (ductile-ish, tail barely matters) to severe (brittle ceramic).
# Sea ice in tension is usually placed at the low end of the brittle range.
M_GRID = [1, 2, 3, 5, 8, 12, 16, 20, 30, 50]
M_REF = 8.0            # representative modulus used for the headline ranking


def scf_eff(scf, vol, m):
    """Volume-weighted m-norm. Computed in log space: scf^50 overflows float64
    for the base case, whose peak SCF is ~20."""
    w = vol / vol.sum()
    a = np.log(scf.clip(min=1e-12))
    peak = a.max()
    # (sum w exp(m a))^(1/m) == exp(peak) * (sum w exp(m(a-peak)))^(1/m)
    return float(np.exp(peak) * np.exp(np.log((w * np.exp(m * (a - peak))).sum()) / m))


def equivalent_percentile(scf, vol, target):
    """The percentile of the (volume-weighted) SCF distribution equal to target."""
    o = np.argsort(scf)
    s, w = scf[o], vol[o]
    c = np.cumsum(w) / w.sum()
    i = np.searchsorted(s, target)
    if i <= 0:
        return 0.0
    if i >= len(s):
        return 100.0
    return float(100.0 * c[i - 1])


def load(dump_dir):
    cases = {}
    for p in sorted(glob.glob(os.path.join(dump_dir, '*.npz'))):
        z = np.load(p, allow_pickle=True)
        name = os.path.splitext(os.path.basename(p))[0]
        mm = re.match(r'WBL_([A-Z]+)_s(\d+)', name)
        tag = mm.group(1) if mm else name
        cases.setdefault(tag, []).append((name, z['scf'], z['vol']))
    return cases


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    dump_dir = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else 'results_weibull'

    cases = load(dump_dir)
    if not cases:
        sys.exit('no .npz dumps found in %s' % dump_dir)

    import csv
    rows = []
    print('%-6s %3s %8s %8s %8s %s' % ('case', 'n', 'P99', 'max', 'm=%g' % M_REF,
                                       '  SCF_eff(m) ->'))
    for tag, reps in sorted(cases.items()):
        per_m = {m: [] for m in M_GRID}
        p99s, maxs, pct_of_ref = [], [], []
        for _, scf, vol in reps:
            for m in M_GRID:
                per_m[m].append(scf_eff(scf, vol, m))
            p99s.append(float(np.percentile(scf, 99)))
            maxs.append(float(scf.max()))
            pct_of_ref.append(equivalent_percentile(scf, vol, scf_eff(scf, vol, M_REF)))
        row = {'case': tag, 'n_packings': len(reps),
               'P99_mean': np.mean(p99s), 'P99_sd': np.std(p99s),
               'max_mean': np.mean(maxs), 'max_sd': np.std(maxs),
               'pct_equiv_at_m%g' % M_REF: np.mean(pct_of_ref)}
        for m in M_GRID:
            row['SCFeff_m%g' % m] = np.mean(per_m[m])
            row['SCFeff_m%g_sd' % m] = np.std(per_m[m])
        rows.append(row)
        print('%-6s %3d %8.3f %8.3f %8.3f   %s' % (
            tag, len(reps), row['P99_mean'], row['max_mean'],
            row['SCFeff_m%g' % M_REF],
            ' '.join('%.2f' % row['SCFeff_m%g' % m] for m in M_GRID)))

    with open(prefix + '.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print('\nwrote %s.csv' % prefix)

    # Does the ranking of cases depend on m? That is the question Section 3.3
    # answers with the P50/P90/P99/max columns; here it is answered continuously.
    print('\nranking by SCF_eff, most concentrated first:')
    for m in M_GRID:
        order = sorted(rows, key=lambda r: -r['SCFeff_m%g' % m])
        print('  m=%-3g %s' % (m, ' > '.join(r['case'] for r in order)))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import os as _os
        import sys as _sys
        _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
        import figstyle as fs
        fs.apply()
    except Exception:
        print('(matplotlib/figstyle unavailable, skipping figure)')
        return

    # Okabe-Ito, matching the rest of the figures; BASE deliberately the one
    # dark saturated colour since it is the case every conclusion turns on.
    COL = {'BASE': fs.VERM, 'CHAN': fs.BLUE, 'GAS': fs.ORANGE,
           'ELON': fs.GREEN, 'POCK': fs.PURPLE, 'CTRL': fs.SKY}
    order = ['BASE', 'CHAN', 'GAS', 'ELON', 'POCK', 'CTRL']
    rows.sort(key=lambda r: order.index(r['case']) if r['case'] in order else 99)

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for r in rows:
        y = [r['SCFeff_m%g' % m] for m in M_GRID]
        e = [r['SCFeff_m%g_sd' % m] for m in M_GRID]
        ax[0].errorbar(M_GRID, y, yerr=e, marker='o', ms=4, capsize=2,
                       color=COL.get(r['case']), label=r['case'])
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Weibull modulus $m$')
    ax[0].set_ylabel(r'$\mathrm{SCF}_{\mathrm{eff}}(m)$')
    ax[0].set_title('effective concentration vs weakest-link severity')
    ax[0].legend(fontsize=7)

    # Second panel: the one ordering that actually depends on m. Ratios to a
    # common reference would obscure it, since every case rises with m together.
    d = {r['case']: r for r in rows}
    if 'CHAN' in d and 'GAS' in d:
        ratio = [d['CHAN']['SCFeff_m%g' % m] / d['GAS']['SCFeff_m%g' % m]
                 for m in M_GRID]
        ax[1].plot(M_GRID, ratio, marker='o', ms=4, color=fs.BLUE,
                   label='channels / gas voids')
        ax[1].axhline(1.0, lw=1.0, color=fs.BLACK, ls='--')
        cross = None
        for i in range(len(M_GRID) - 1):
            if (ratio[i] - 1.0) * (ratio[i + 1] - 1.0) < 0:
                cross = M_GRID[i] + (M_GRID[i + 1] - M_GRID[i]) * \
                    (1.0 - ratio[i]) / (ratio[i + 1] - ratio[i])
        if cross:
            ax[1].axvline(cross, lw=1.0, color=fs.VERM, ls=':')
            ax[1].annotate('order reverses at $m\\approx%.0f$' % cross,
                           xy=(cross, 1.0), xytext=(cross * 1.3, 0.92),
                           color=fs.VERM, fontsize=11.5,
                           arrowprops=dict(arrowstyle='->', color=fs.VERM, lw=1.2))
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Weibull modulus $m$')
    ax[1].set_ylabel(r'$\mathrm{SCF}_{\mathrm{eff}}$ ratio')
    ax[1].set_title(r'the one ordering that depends on $m$')
    ax[1].legend(loc='upper left', fontsize=11)
    fig.tight_layout()
    fig.savefig(prefix + '.png', dpi=200)
    print('wrote %s.png' % prefix)


if __name__ == '__main__':
    main()
