r"""Generate the figure set for ice_rve.tex, from the campaign results.

One script so the figures stay consistent with each other and reproducible
from the committed CSVs. Schematics are not here: Figure 1 is ported from
main_fix and stays in tikz.

Each figure carries one claim, and where a prediction was made before the
campaign ran, the prediction is DRAWN rather than fitted so the reader can see
it was not tuned afterwards.

    python3 analysis/make_ice_rve_figs.py [outdir]
"""
import csv
import os
import re
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import layered_law as law

RES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'results')


def load(name, pattern, cols=('E_x',)):
    """Group a results file by run-id fields matched from `pattern`."""
    out = defaultdict(list)
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        return None
    for r in csv.DictReader(open(path, encoding='utf8', errors='replace')):
        m = re.match(pattern, r.get('run_id', ''))
        if not m:
            continue
        try:
            vals = tuple(float(r[c]) / 1e9 for c in cols)
        except (ValueError, KeyError, TypeError):
            continue
        out[m.groups()].append(vals)
    return out


def agg(g, key_fn):
    """Mean and population s.d. per key, sorted."""
    d = defaultdict(list)
    for k, v in g.items():
        d[key_fn(k)].extend(v)
    ks = sorted(d)
    arr = [np.array(d[k]) for k in ks]
    return (np.array(ks, dtype=float),
            np.array([a[:, 0].mean() for a in arr]),
            np.array([a[:, 0].std() for a in arr]))


# ---------------------------------------------------------------- Fig 3
def fig_rve(outdir):
    """The layered cell is a representative volume -- once the microstructure
    is actually held fixed."""
    dens = load('results_bracket_density.csv', r'BRKD_L(\d+)_(und|drn)_s\d')
    spac = load('results_bracket_spacing.csv', r'BRKS_L(\d+)_(und|drn)_s\d')
    if not dens or not spac:
        print('  fig3: missing data, skipped')
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.2), sharey=True)
    for a, g, letter, note in (
            (ax[0], dens, 'a',
             'bridge density held\nconverged: CV 1.9% / 0.8%'),
            (ax[1], spac, 'b',
             'bridge count held instead' '\n'
             r'drifts as $L^{-1.14}$, CV 39%')):
        for state, c, mk, lab in (('drn', fs.BLUE, 'o', 'drained'),
                                  ('und', fs.VERM, 's', 'undrained')):
            sub = {k: v for k, v in g.items() if k[1] == state}
            L, m, s = agg(sub, lambda k: int(k[0]) / 1000.0)
            a.errorbar(L, m, yerr=s, marker=mk, color=c, capsize=3, label=lab)
        a.set_xlabel(r'cell edge $L$   [model units]')
        a.text(0.04, 0.06, note, transform=a.transAxes, fontsize=11,
               color='0.25')
        a.text(0.015, 0.965, '(%s)' % letter, transform=a.transAxes,
               fontsize=13, fontweight='bold', va='top')
        a.set_yscale('log'); a.set_ylim(0.15, 9)
        a.legend(fontsize=10.5, loc='center left')
    ax[0].set_ylabel(r'$E_x$   [GPa]')
    fig.tight_layout()
    p = os.path.join(outdir, 'fig3_rve.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)


# ---------------------------------------------------------------- Fig 5
def fig_elimination(outdir):
    """Five candidate explanations, each measured, against the gap."""
    # every value here is measured and its provenance is in the commit log
    rows = [('inclusion aspect ratio\nneedle to 4:1 plate', 1.03),
            ('brine channels\nadded at matched fraction', 1.00),
            ('pore-pressure drainage\npocket morphology', 1.04),
            ('inversion artefact\nlinear-profile assumption', 1.057),
            ('porosity route\n(needs $\\phi=0.52$ vs 0.23 available)', 1.00)]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    y = np.arange(len(rows))
    ax.barh(y, [r[1] for r in rows], height=0.55, color=fs.SKY,
            edgecolor=fs.BLUE)
    ax.axvspan(4.0, 7.0, color=fs.VERM, alpha=0.22, zorder=0)
    ax.text(5.5, len(rows) - 0.35, 'factor needing\nexplanation', fontsize=11.5,
            color=fs.VERM, ha='center', va='top')
    for i, (_, v) in enumerate(rows):
        ax.text(v * 1.03, i, r'$\times$%.2f' % v, va='center', fontsize=11)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=10.5)
    ax.set_xscale('log'); ax.set_xlim(0.9, 9)
    ax.set_xlabel('factor by which the transverse modulus changes')
    ax.grid(axis='y', alpha=0)
    fig.tight_layout()
    p = os.path.join(outdir, 'fig5_elimination.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)


# ---------------------------------------------------------------- Fig 6
def fig_confinement(outdir):
    """Releasing the fill's bulk modulus: decisive for layers, not for pockets,
    and invariant along the layers, which is what names the mechanism."""
    K = np.array([2.2, 0.22, 0.022, 0.0022])
    Ex = np.array([4.605, 1.612, 0.780, 0.671])
    Ez = np.array([7.493, 7.487, 7.476, 7.520])
    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.2))
    a = ax[0]
    a.plot(K, Ex, 'o-', color=fs.BLUE, ms=9, label=r'$E_x$, across the layers')
    a.plot(K, Ez, 's--', color=fs.GREEN, ms=9, label=r'$E_z$, along the layers')
    a.annotate(r'$\times%.1f$' % (Ex[0] / Ex[-1]), xy=(0.0022, 0.9),
               xytext=(0.02, 0.30), fontsize=12, color=fs.BLUE,
               arrowprops=dict(arrowstyle='<->', color=fs.BLUE, lw=1.4))
    a.text(0.004, 8.6, r'$E_z$ flat to 0.6%', fontsize=11.5, color=fs.GREEN)
    a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(0.4, 14)
    a.set_xlabel(r"fill bulk modulus $K$   [GPa]")
    a.set_ylabel(r'$E$   [GPa]')
    a.text(0.015, 0.965, '(a)', transform=a.transAxes,
            fontsize=13, fontweight='bold', va='top')
    a.legend(fontsize=10.5, loc='lower right')

    b = ax[1]
    names = ['pockets', 'layers']
    ratio = [6.110 / 5.903, Ex[0] / Ex[-1]]
    b.bar(names, ratio, width=0.5, color=[fs.ORANGE, fs.BLUE], alpha=0.9)
    for i, v in enumerate(ratio):
        b.text(i, v * 1.04, r'$\times$%.2f' % v, ha='center', fontsize=12)
    b.set_yscale('log'); b.set_ylim(0.9, 12)
    b.set_ylabel('undrained / drained')
    b.text(0.015, 0.965, '(b)', transform=b.transAxes,
            fontsize=13, fontweight='bold', va='top')
    b.grid(axis='x', alpha=0)
    fig.tight_layout()
    p = os.path.join(outdir, 'fig6_confinement.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)


# ---------------------------------------------------------------- Fig 7
def fig_constriction(outdir):
    """Subdividing a fixed bridge area. The N^0.5 line is the prediction,
    drawn and not fitted."""
    g = load('results_bracket_nbridges.csv', r'BRKG_n(\d+)_(und|drn)_s\d')
    if not g:
        print('  fig7: missing data, skipped')
        return
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    for state, c, mk, lab in (('drn', fs.BLUE, 'o', 'drained'),
                              ('und', fs.VERM, 's', 'undrained')):
        sub = {k: v for k, v in g.items() if k[1] == state}
        N, m, s = agg(sub, lambda k: int(k[0]))
        ax.errorbar(N, m, yerr=s, marker=mk, color=c, ms=9, capsize=3,
                    label=lab)
        p = np.polyfit(np.log(N), np.log(m), 1)[0]
        ax.text(N[-1] * 1.12, m[-1], r'$N^{%.3f}$' % p, color=c, fontsize=12,
                va='center')
    sub = {k: v for k, v in g.items() if k[1] == 'drn'}
    N, m, _ = agg(sub, lambda k: int(k[0]))
    ax.plot(N, m[0] * (N / N[0]) ** 0.5, ':', color=fs.BLACK, lw=2.0,
            label=r'prediction $N^{1/2}$ (spreading compliance)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('bridges per layer, at fixed total bridge area')
    ax.set_ylabel(r'$E_x$   [GPa]')
    ax.legend(fontsize=10.5, loc='center left')
    fig.tight_layout()
    p = os.path.join(outdir, 'fig7_constriction.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)


# ---------------------------------------------------------------- Fig 8
def fig_spacing(outdir):
    """Modulus against lamellar spacing, with the measured spacing marked."""
    g = load('results_bracket_nlayers.csv', r'BRKN_n(\d+)_(und|drn)_s\d')
    if not g:
        print('  fig8: missing data, skipped')
        return
    L_CELL, CELL_MM = 0.5, 3.0
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for state, c, mk, lab in (('drn', fs.BLUE, 'o', 'drained'),
                              ('und', fs.VERM, 's', 'undrained')):
        sub = {k: v for k, v in g.items() if k[1] == state}
        n, m, s = agg(sub, lambda k: int(k[0]))
        a0_mm = L_CELL / n * (CELL_MM / L_CELL)
        ax.errorbar(a0_mm, m, yerr=s, marker=mk, color=c, ms=9, capsize=3,
                    label=lab)
        p = np.polyfit(np.log(a0_mm), np.log(m), 1)[0]
        ax.text(a0_mm[0] * 1.06, m[0], r'$a_0^{%+.2f}$' % p, color=c,
                fontsize=12, va='center')
    ax.axvspan(0.20, 0.50, color=fs.GREEN, alpha=0.20, zorder=0)
    ax.text(0.31, 0.62, 'measured\n(Pringle)', fontsize=11, color=fs.GREEN,
            ha='center')
    ax.annotate('solved cells stop here', xy=(0.75, 0.42), xytext=(1.1, 0.30),
                fontsize=10.5, color='0.3',
                arrowprops=dict(arrowstyle='->', color='0.4', lw=1.2))
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'lamellar spacing $a_0$   [mm]')
    ax.set_ylabel(r'$E_x$   [GPa]')
    ax.legend(fontsize=10.5, loc='upper left')
    fig.tight_layout()
    p = os.path.join(outdir, 'fig8_spacing.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)


# ---------------------------------------------------------------- Fig 9
def fig_bridge(outdir):
    """Bridge fraction at fixed brine content. The drained cell follows a power
    law in b; the undrained one does not respond to b at all. Both had to be
    measured, because the confinement argument predicts the second and the
    Assur plane-of-weakness argument predicts only the first."""
    g = load('results_bracket_bridge.csv', r'BRKB_b(\d+)_(und|drn)_s\d')
    if not g:
        print('  fig9: missing data, skipped')
        return
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.3))
    a = ax[0]
    fits = {}
    for state, c, mk, lab in (('drn', fs.BLUE, 'o', 'drained'),
                              ('und', fs.VERM, 's', 'undrained')):
        sub = {k: v for k, v in g.items() if k[1] == state}
        b, m, s = agg(sub, lambda k: int(k[0]) / 1000.0)
        a.errorbar(b, m, yerr=s, marker=mk, color=c, ms=9, capsize=3, label=lab)
        p, cov = np.polyfit(np.log(b), np.log(m), 1, cov=True)
        fits[state] = (p[0], np.sqrt(cov[0, 0]))
        a.text(b[-1] * 1.1, m[-1], r'$b^{%.2f}$' % p[0], color=c, fontsize=12,
               va='center')
    a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(0.15, 12)
    a.set_xlabel(r'ice fraction of the layer plane $b$')
    a.set_ylabel(r'$E_x$   [GPa]')
    a.text(0.015, 0.965, '(a)', transform=a.transAxes, fontsize=13,
           fontweight='bold', va='top')
    a.legend(fontsize=10.5, loc='upper left')

    # the drainage ratio: largest exactly where the bridges are sparsest
    b2, md, _ = agg({k: v for k, v in g.items() if k[1] == 'drn'},
                    lambda k: int(k[0]) / 1000.0)
    _, mu, _ = agg({k: v for k, v in g.items() if k[1] == 'und'},
                   lambda k: int(k[0]) / 1000.0)
    r = mu / md
    bb = ax[1]
    bb.plot(b2, r, 'D-', color=fs.PURPLE, ms=9)
    for x, y in zip(b2, r):
        bb.text(x, y * 1.10, r'$\times%.1f$' % y, ha='center', fontsize=10.5)
    bb.set_xscale('log'); bb.set_yscale('log'); bb.set_ylim(2, 32)
    bb.set_xlabel(r'ice fraction of the layer plane $b$')
    bb.set_ylabel('undrained / drained')
    bb.text(0.015, 0.965, '(b)', transform=bb.transAxes, fontsize=13,
            fontweight='bold', va='top')
    bb.text(0.55, 0.80, 'drainage matters most\nwhere the bridges are sparsest',
            transform=bb.transAxes, fontsize=10.5, color='0.25', ha='center')
    fig.tight_layout()
    p = os.path.join(outdir, 'fig9_bridge.png')
    fig.savefig(p, dpi=170); print('  wrote %s' % p)
    print('     drained  n = %.3f +- %.3f' % fits['drn'])
    print('     undrained n = %.3f +- %.3f' % fits['und'])


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else RES
    print('writing ice_rve figures into %s' % outdir)
    fig_rve(outdir)
    fig_elimination(outdir)
    fig_confinement(outdir)
    fig_constriction(outdir)
    fig_spacing(outdir)
    fig_bridge(outdir)
    print('Fig 1 is ported from main_fix (tikz); Fig 2 needs mesh renders;')
    print('Figs 9 and 10 come from plot_ez_closure.py and plot_match_ez.py.')


if __name__ == '__main__':
    main()
