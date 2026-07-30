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
        import figstyle  # noqa: F401  (shared rc)
    except Exception:
        print('(matplotlib/figstyle unavailable, skipping figure)')
        return

    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.0))
    for r in rows:
        y = [r['SCFeff_m%g' % m] for m in M_GRID]
        e = [r['SCFeff_m%g_sd' % m] for m in M_GRID]
        ax[0].errorbar(M_GRID, y, yerr=e, marker='o', ms=3, capsize=2, label=r['case'])
    ax[0].set_xscale('log')
    ax[0].set_xlabel('Weibull modulus $m$')
    ax[0].set_ylabel(r'$\mathrm{SCF}_{\mathrm{eff}}(m)$')
    ax[0].set_title('effective concentration vs weakest-link severity')
    ax[0].legend(fontsize=7)

    base = [r for r in rows if r['case'] == 'BASE']
    ref = base[0] if base else rows[-1]
    for r in rows:
        y = [r['SCFeff_m%g' % m] / ref['SCFeff_m%g' % m] for m in M_GRID]
        ax[1].plot(M_GRID, y, marker='o', ms=3, label=r['case'])
    ax[1].set_xscale('log')
    ax[1].axhline(1.0, lw=0.8, color='0.5')
    ax[1].set_xlabel('Weibull modulus $m$')
    ax[1].set_ylabel(r'$\mathrm{SCF}_{\mathrm{eff}}$ / %s' % ref['case'])
    ax[1].set_title('ranking stability')
    fig.tight_layout()
    fig.savefig(prefix + '.png', dpi=200)
    print('wrote %s.png' % prefix)


if __name__ == '__main__':
    main()
