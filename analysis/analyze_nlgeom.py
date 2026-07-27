"""Study #8 analysis: large-deformation (nlgeom) homogenization.

Reads the reaction-based nominal stress-strain paths (curves_nlgeom_{lin,ten,cmp}.csv)
and the moduli summaries (results_nlgeom_*.csv) for three column slices, each loaded
on ONE shared mesh three ways: linear reference (small-strain theory), finite-strain
tension, finite-strain compression. Builds study_nlgeom.png (matplotlib only).

Findings:
  * the linear reference is exactly straight (sec/E0 = 1.0000) -- validates the
    RP-reaction extraction;
  * COMPRESSION converges to 2% for every slice and stiffens uniformly (+1.5-1.8%);
  * TENSION softens (z25: -1.7% at 2%) but the percolated channelled base is
    geometrically UNSTABLE in tension -- the solve loses convergence early
    (z65 at 0.75%, z95 at 0.1% strain), a tension/compression asymmetry in both
    modulus and stability that is absent from the small-strain homogenization.
"""
import csv
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Shared enlarged fonts + Okabe-Ito palette (see figstyle.py). This figure is
# stress-strain, not depth-resolved, so it keeps a linear x-axis (no z/H flip).
from figstyle import BLUE, VERM
from figstyle import apply as _apply_style
_apply_style()

SLICES = ['z25', 'z65', 'z95']
TITLES = {'z25': '(a) z/H=0.25  low porosity',
          'z65': '(b) z/H=0.65  transition',
          'z95': '(c) z/H=0.95  channelled base'}

def load_curves(case):
    d = defaultdict(list)
    for r in csv.DictReader(open('curves_nlgeom_%s.csv' % case)):
        d[r['run_id']].append((float(r['eps']), float(r['sigma'])))
    for k in d:
        d[k].sort()
    return d

def load_summary(case):
    return {r['run_id']: r for r in csv.DictReader(open('results_nlgeom_%s.csv' % case))}

def main():
    lin, ten, cmp = load_curves('lin'), load_curves('ten'), load_curves('cmp')
    sl, st, sc = load_summary('lin'), load_summary('ten'), load_summary('cmp')

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.4))
    for ax, s in zip(axes, SLICES):
        E0 = float(sl['NLGLIN_%s' % s]['E0'])
        # linear reference line through origin, slope E0 (both branches)
        xr = [-0.021, 0.021]
        ax.plot([x * 100 for x in xr], [E0 * x / 1e6 for x in xr],
                '--', color='0.55', lw=1.6, label='linear (small-strain)', zorder=1)
        # compression branch (eps < 0)
        cc = cmp['NLGCMP_%s' % s]
        ax.plot([e * 100 for e, _ in cc], [sig / 1e6 for _, sig in cc],
                'o-', color=BLUE, ms=4, label='nlgeom compression', zorder=3)
        # tension branch (eps > 0), mark the last converged point if it stalled
        tt = ten['NLGTEN_%s' % s]
        ax.plot([e * 100 for e, _ in tt], [sig / 1e6 for _, sig in tt],
                'o-', color=VERM, ms=4, label='nlgeom tension', zorder=3)
        eps_reached = tt[-1][0]
        if eps_reached < 0.019:                     # tension lost stability early
            ax.plot(eps_reached * 100, tt[-1][1] / 1e6, 'x', color='k', ms=12,
                    mew=2.8, zorder=4)
            ax.annotate('tensile\ninstability\n(%.1f%% strain)' % (eps_reached * 100),
                        xy=(eps_reached * 100, tt[-1][1] / 1e6),
                        xytext=(0.30, 0.30), textcoords='axes fraction',
                        fontsize=11, ha='center',
                        arrowprops=dict(arrowstyle='->', color='k', lw=1.2))
        # annotate the geometric correction where each branch converged
        secC = float(sc['NLGCMP_%s' % s]['sec_over_E0'])
        txt = 'compr: %+.1f%%' % (100 * (secC - 1))
        if eps_reached >= 0.019:
            secT = float(st['NLGTEN_%s' % s]['sec_over_E0'])
            txt = 'tens: %+.1f%%   ' % (100 * (float(st['NLGTEN_%s' % s]['sec_over_E0']) - 1)) + txt
        ax.text(0.5, 0.02, txt, transform=ax.transAxes, ha='center', va='bottom',
                fontsize=11, bbox=dict(boxstyle='round', fc='w', ec='0.7', alpha=0.9))
        ax.axhline(0, color='0.8', lw=0.8); ax.axvline(0, color='0.8', lw=0.8)
        ax.set_title(TITLES[s])
        ax.set_xlabel('Nominal strain $\\varepsilon_{nom}$ (%)')
        ax.set_ylabel('Nominal stress $\\sigma_{nom}$ (MPa)')
        if s == 'z25':
            ax.legend(fontsize=10.5, loc='upper left')
    fig.tight_layout()
    fig.savefig('study_nlgeom.png', dpi=200)
    print('wrote study_nlgeom.png')

    print('\n%-6s %10s %12s %12s %14s' % ('slice', 'E0(GPa)', 'tens sec/E0',
                                          'compr sec/E0', 'tens eps_max'))
    for s in SLICES:
        E0 = float(sl['NLGLIN_%s' % s]['E0']) / 1e9
        secT = float(st['NLGTEN_%s' % s]['sec_over_E0'])
        secC = float(sc['NLGCMP_%s' % s]['sec_over_E0'])
        emx = float(st['NLGTEN_%s' % s]['eps_max'])
        print('%-6s %10.3f %12.4f %12.4f %14.4f' % (s, E0, secT, secC, emx))

if __name__ == '__main__':
    main()
