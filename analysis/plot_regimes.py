r"""The four regimes of the closure, drawn.

Columnar sea ice grows as vertical ice platelets separated by brine layers, and
the plate spacing a_0 is fixed at the growth interface. So the lamellar planes
exist at every depth: nothing about the geometry switches on as the ice warms.
What changes is how much brine sits in those planes and, crucially, how it is
connected. Each of Pringle's three percolation thresholds marks one connectivity
change, and each connectivity change alters the load path.

Layout: a 2x2 grid, one regime per panel. Within each panel, the left half is a
vertical section through three platelets and the two brine planes between them,
the right half is one of those planes seen face on (the view b measures). Below
the grid, the closure itself, shaded by the same four regimes.

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
     r'$E=E_{\rm ice}(1-1.65\,\phi)$'),
    (r'$0.046 < \phi < 0.09$', 'drained pockets',
     r'same, $\div 1.04$'),
    (r'$0.09 < \phi < 0.14$', 'bridge constrictions',
     r'$\times\, b^{\,n(b)\,w}$'),
    (r'$\phi > 0.14$', 'strut network',
     r'same law (measured $b^{\,0.84-0.88}$)'),
]


def platelets(ax, mode):
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
                    arrowprops=dict(arrowstyle='->', color=fs.VERM, lw=1.8))
        ax.text(0.42, 0.995, 'drains', fontsize=11, color=fs.VERM,
                va='top')
    if mode == 'struts':
        ax.annotate('', xy=(0.97, 0.48), xytext=(0.03, 0.48),
                    arrowprops=dict(arrowstyle='<->', color=fs.VERM, lw=1.7))
        ax.text(0.5, 0.55, 'brine crosses', fontsize=11, color=fs.VERM,
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
        lab = r'$b<0.16$: sparse bridges'
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(lab, fontsize=11.5, labelpad=4)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    fig = plt.figure(figsize=(12.6, 10.6))

    # Four regimes as a 2x2 grid (each panel roughly twice the area of the old
    # 1x4 strip), with the closure below spanning the full width.
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.28],
                          hspace=0.42, wspace=0.14)
    modes = ['sealed', 'drained', 'bridge', 'struts']
    for i, (rng_lab, name, law) in enumerate(REGIMES):
        row, col = divmod(i, 2)
        cell = gs[row, col].subgridspec(2, 2, height_ratios=[0.30, 1.0],
                                        wspace=0.08, hspace=0.10)
        axt = fig.add_subplot(cell[0, :])
        axt.axis('off')
        axt.text(0.0, 0.74, r'%s\ \ \ %s' % (rng_lab, name),
                 fontsize=14.5, fontweight='bold', va='center')
        axt.text(0.0, 0.14, law, fontsize=13.5, color=fs.BLUE, va='center')

        axv = fig.add_subplot(cell[1, 0])
        platelets(axv, modes[i])
        axp = fig.add_subplot(cell[1, 1])
        inplane(axp, modes[i])

    # The closure itself, shaded by the same four regimes as the pictures
    # above, so the reader meets the mechanism and the law on one page.
    axc = fig.add_subplot(gs[2, :])
    phi = np.linspace(0.002, 0.235, 3000)
    lo, mid, hi = ez.E_band(phi, floor=0.0)
    edges = [0.0, ez.PHI_DRAIN, ez.PHI_LAYER, ez.PHI_CROSS, ez.PHI_0, 0.235]
    shades = ['#f5f5f5', '#e8eef4', '#dbe6f0', '#c9d9ea', '#efe3e3']
    for a, b, c in zip(edges[:-1], edges[1:], shades):
        axc.axvspan(a, b, color=c, zorder=0)
    axc.axvspan(ez.PHI_DRAIN - ez.PHI_DRAIN_SD, ez.PHI_DRAIN + ez.PHI_DRAIN_SD,
                color=fs.VERM, alpha=0.20, zorder=0)
    axc.fill_between(phi, lo, hi, color=fs.SKY, alpha=0.40,
                     label=r'closure, $n(b)$, $\pm2\,\mathrm{rms}$')
    axc.plot(phi, np.maximum(mid, 0.04), color=fs.BLUE, lw=2.8,
             label=r'closure, $n(b)$')
    axc.plot(phi, ez.E_ICE * (1 - 1.65 * phi), color=fs.ORANGE, lw=2.0,
             ls='--', label='pocket cells (measured)')
    axc.plot(phi, 9.5 * (1 - np.sqrt(phi)) ** 4, color=fs.BLACK, lw=1.6,
             label='Weeks & Assur 1967')
    for x, lab in ((ez.PHI_DRAIN, r'$0.046$'), (ez.PHI_LAYER, r'$0.09$'),
                   (ez.PHI_CROSS, r'$0.14$'), (ez.PHI_0, r'$\phi_0=0.20$')):
        axc.axvline(x, color='0.4', lw=1.0, ls=':')
        axc.text(x, 17.0, lab, fontsize=12, color='0.3', ha='center', va='top')
    axc.text(0.216, 0.30, 'skeletal:\nno model here', fontsize=11.5,
             color=fs.VERM, ha='center')
    axc.set_yscale('log'); axc.set_ylim(0.04, 22); axc.set_xlim(0, 0.235)
    axc.set_xlabel(r'brine volume fraction $\phi$', fontsize=15)
    axc.set_ylabel(r'$E$   [GPa]', fontsize=15)
    axc.legend(loc='lower left', fontsize=12.5)

    fig.subplots_adjust(left=0.055, right=0.985, top=0.96, bottom=0.055)
    p = os.path.join(out, 'fig_regimes.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
