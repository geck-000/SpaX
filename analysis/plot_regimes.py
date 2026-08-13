r"""The four regimes of the closure, drawn.

Columnar sea ice grows as vertical ice platelets separated by brine layers, and
the plate spacing a_0 is fixed at the growth interface. So the lamellar planes
exist at every depth: nothing about the geometry switches on as the ice warms.
What changes is how much brine sits in those planes and, crucially, how it is
connected. Each of Pringle's three percolation thresholds marks one connectivity
change, and each connectivity change alters the load path.

Top row: a vertical section through three platelets and the two brine planes
between them, plus the in-plane view of one plane, which is what b measures.
Bottom: the closure, with the regimes shaded.

    python3 plot_regimes.py [outdir]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse, Circle

import ez_closure as ez

ICE = '#dce9f4'
ICE_E = '#8fb4d0'
BRINE = '#E69F00'

REGIMES = [
    (r'$\phi < 0.046$', 'sealed pockets',
     'brine isolated in every direction;\nice continuous; pore pressure sealed,\n'
     'so brine resists at its bulk modulus',
     r'$E=E_{\rm ice}(1-1.65\phi)$'),
    (r'$0.046 < \phi < 0.09$', 'drained pockets',
     'brine now connects VERTICALLY and\ncan escape; pockets still separate\n'
     'within a plane, so ice still continuous',
     r'same, $\div 1.04$'),
    (r'$0.09 < \phi < 0.14$', 'bridge constrictions',
     'brine spans each PLANE; load crossing\nit must funnel through ice bridges,\n'
     'so stress spreads into each contact',
     r'$\times\, b^{\,n_{\rm eff}\,w}$, $n\!\approx\!0.5$'),
    (r'$\phi > 0.14$', 'strut bending',
     'brine crosses BETWEEN planes too:\nplatelets breached, ice left as sparse\n'
     'struts that carry load by bending',
     r'$\times\, b^{2}$'),
]


def platelets(ax, brine_frac, mode):
    """Vertical section: three ice platelets, two brine planes between them."""
    ax.add_patch(Rectangle((0, 0), 1, 1, fc='white', ec='none'))
    xs = [0.0, 0.36, 0.64, 1.0]          # platelet edges
    for x0, x1 in ((0.0, 0.30), (0.36, 0.64), (0.70, 1.0)):
        ax.add_patch(Rectangle((x0, 0.06), x1 - x0, 0.88, fc=ICE,
                               ec=ICE_E, lw=0.9))
    for xc in (0.33, 0.67):
        w = 0.055
        if mode in ('sealed', 'drained'):
            # brine necked into discrete pockets inside the plane
            for yc in (0.20, 0.42, 0.64, 0.84):
                ax.add_patch(Ellipse((xc, yc), w * 1.5, 0.11, fc=BRINE,
                                     ec='none'))
        elif mode == 'bridge':
            # continuous brine sheet, interrupted by a few ice bridges
            ax.add_patch(Rectangle((xc - w / 2, 0.06), w, 0.88, fc=BRINE,
                                   ec='none'))
            for yc in (0.30, 0.70):
                ax.add_patch(Rectangle((xc - w / 2, yc - 0.055), w, 0.11,
                                       fc=ICE, ec=ICE_E, lw=0.7))
        else:  # struts
            ax.add_patch(Rectangle((xc - w / 2, 0.06), w, 0.88, fc=BRINE,
                                   ec='none'))
            ax.add_patch(Rectangle((xc - w / 2, 0.48), w, 0.045, fc=ICE,
                                   ec=ICE_E, lw=0.7))
            # brine has also broken through the platelets themselves
            for (px, py) in ((0.17, 0.30), (0.17, 0.66), (0.50, 0.48),
                             (0.85, 0.34), (0.85, 0.70)):
                ax.add_patch(Ellipse((px, py), 0.10, 0.07, fc=BRINE,
                                     ec='none', alpha=0.95))
    if mode == 'drained':
        ax.annotate('', xy=(0.33, 0.99), xytext=(0.33, 0.02),
                    arrowprops=dict(arrowstyle='->', color=fs.VERM, lw=1.6))
        ax.text(0.40, 0.995, 'drains', fontsize=7.5, color=fs.VERM,
                va='top')
    if mode == 'struts':
        ax.annotate('', xy=(0.97, 0.48), xytext=(0.03, 0.48),
                    arrowprops=dict(arrowstyle='<->', color=fs.VERM, lw=1.5))
        ax.text(0.5, 0.53, 'brine crosses', fontsize=7.5, color=fs.VERM,
                ha='center')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.06)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def inplane(ax, mode):
    """Face-on view of ONE brine plane. Ice fraction here is b."""
    ax.add_patch(Rectangle((0, 0), 1, 1, fc=BRINE, ec='0.5', lw=0.8))
    rng = np.random.default_rng(4)
    if mode in ('sealed', 'drained'):
        ax.add_patch(Rectangle((0, 0), 1, 1, fc=ICE, ec='0.5', lw=0.8))
        for _ in range(9):
            ax.add_patch(Ellipse(rng.uniform(0.12, 0.88, 2), 0.13, 0.09,
                                 fc=BRINE, ec='none'))
        lab = (r'$b\to1$: ice spans the plane'
               if mode == 'sealed' else
               r'$b\to1$ still: drainage is a'
               '\n' r'vertical change, not an in-plane one')
    elif mode == 'bridge':
        # non-overlapping bridges on a coarse jittered lattice
        pts = [(0.25, 0.25), (0.72, 0.24), (0.24, 0.72), (0.70, 0.71),
               (0.48, 0.48)]
        for px, py in pts:
            ax.add_patch(Circle((px + rng.uniform(-0.03, 0.03),
                                 py + rng.uniform(-0.03, 0.03)), 0.125,
                                fc=ICE, ec=ICE_E, lw=0.8))
        lab = r'$b\approx0.16$--$0.33$: bridges'
    else:
        pts = [(0.27, 0.28), (0.71, 0.26), (0.26, 0.71), (0.72, 0.73),
               (0.49, 0.50)]
        for px, py in pts:
            ax.add_patch(Circle((px, py), 0.052, fc=ICE, ec=ICE_E, lw=0.8))
        lab = r'$b<0.16$: sparse struts'
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(lab, fontsize=7.4, labelpad=2)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    fig = plt.figure(figsize=(13.6, 9.6))
    gs = fig.add_gridspec(4, 4, height_ratios=[1.35, 0.95, 0.72, 1.45],
                          hspace=0.42, wspace=0.22)
    modes = ['sealed', 'drained', 'bridge', 'struts']
    for i, (rng_lab, name, why, law) in enumerate(REGIMES):
        ax = fig.add_subplot(gs[0, i])
        platelets(ax, None, modes[i])
        ax.set_title('%s\n%s' % (rng_lab, name), fontsize=10.5,
                     fontweight='bold', pad=6)
        ax2 = fig.add_subplot(gs[1, i])
        inplane(ax2, modes[i])
        # description and law get their own row so nothing is clipped
        axt = fig.add_subplot(gs[2, i])
        axt.axis('off')
        axt.text(0.5, 0.98, why, fontsize=8.2, color='0.25', ha='center',
                 va='top', linespacing=1.35)
        axt.text(0.5, 0.06, law, fontsize=9.4, color=fs.BLUE, ha='center',
                 va='bottom')

    ax = fig.add_subplot(gs[3, :])
    phi = np.linspace(1e-4, 0.235, 3000)
    E = ez.E_of_phi(phi, floor=0.0)
    edges = [0.0, ez.PHI_DRAIN, ez.PHI_LAYER, ez.PHI_CROSS, ez.PHI_0, 0.235]
    shades = ['#f5f5f5', '#e8eef4', '#dbe6f0', '#c9d9ea', '#efe3e3']
    for lo, hi, c in zip(edges[:-1], edges[1:], shades):
        ax.axvspan(lo, hi, color=c, zorder=0)
    ax.axvspan(ez.PHI_DRAIN - ez.PHI_DRAIN_SD, ez.PHI_DRAIN + ez.PHI_DRAIN_SD,
               color=fs.VERM, alpha=0.20, zorder=0)
    ax.plot(phi, np.maximum(E, 0.04), color=fs.BLUE, lw=2.8)
    for x, lab in ((ez.PHI_DRAIN, r'$0.046$' '\n' 'vertical'),
                   (ez.PHI_LAYER, r'$0.09$' '\n' 'in-plane'),
                   (ez.PHI_CROSS, r'$0.14$' '\n' 'across-plane'),
                   (ez.PHI_0, r'$0.20$' '\n' r'$\phi_0$')):
        ax.axvline(x, color='0.4', lw=1.0, ls=':')
        ax.text(x, 15.5, lab, fontsize=8.6, color='0.3', ha='center',
                va='top')
    ax.text(0.215, 0.35, 'skeletal:\nno model here', fontsize=8.4,
            color=fs.VERM, ha='center')
    ax.set_yscale('log'); ax.set_ylim(0.04, 22); ax.set_xlim(0, 0.235)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.set_ylabel(r'$E$   [GPa]')
    fig.subplots_adjust(left=0.06, right=0.985, top=0.93, bottom=0.07)
    p = os.path.join(out, 'fig_regimes.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
