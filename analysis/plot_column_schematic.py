r"""The sheet, the depth coordinate, and which description applies where.

Drawn rather than rendered from a mesh, because the point is what stays the
same down the column and what changes, and a photograph of one cell cannot show
that. The lamellae run the full height at a spacing fixed when the ice grew;
only the brine inside them changes with depth.

The two cells beside the column are schematic and correct for their depth --
the earlier version of this figure called out a pocket-and-channel mesh for the
base slice, which is the morphology the base no longer has.

    python3 analysis/plot_column_schematic.py [outdir]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle, Ellipse

import ez_closure as ez

T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])

ICE = '#dcebf5'
ICE_WARM = '#9cc3e4'
BRINE = '#E69F00'
CRIT = '#D55E00'


def thresholds():
    z = (np.arange(10) + 0.5) / 10.0
    phi = ez.brine_volume(T, S)
    zi = np.linspace(0, 1, 4001)
    pi = np.interp(zi, z, phi)
    return (zi[np.argmax(pi >= ez.PHI_DRAIN)],
            zi[np.argmax(pi >= ez.PHI_LAYER)])


def column(ax, z_drain, z_layer):
    """The sheet in section, lamellae at fixed spacing top to bottom."""
    n = 240
    for i in range(n):                       # cold-to-warm wash
        y0, y1 = i / n, (i + 1) / n
        t = i / (n - 1.0)
        c = tuple((1 - t) * np.array(to_rgb(ICE))
                  + t * np.array(to_rgb(ICE_WARM)))
        ax.add_patch(Rectangle((0, y0), 1, y1 - y0, fc=c, ec='none', zorder=0))
    ax.add_patch(Rectangle((0, 0), 1, 1, fc='none', ec='0.25', lw=1.4, zorder=4))

    xs = [0.2, 0.4, 0.6, 0.8]
    for x in xs:
        ax.plot([x, x], [0, 1], color='0.55', lw=0.7, zorder=2)

    # sealed pockets: isolated, above the drainage threshold
    rng = np.random.default_rng(4)
    for x in xs:
        for y in np.arange(0.06, z_drain - 0.01, 0.062):
            ax.add_patch(Ellipse((x, y + rng.uniform(-.008, .008)),
                                 0.030, 0.017, fc=BRINE, ec='none', zorder=3))
    # drained but still separate in the plane
    for x in xs:
        ax.plot([x, x], [z_drain, z_layer], color=BRINE, lw=2.2, zorder=3,
                solid_capstyle='butt', alpha=0.85)
    # spanned plane: continuous brine, ice bridges left in it
    for x in xs:
        ax.plot([x, x], [z_layer, 1.0], color=BRINE, lw=5.0, zorder=3,
                solid_capstyle='butt')
        for y in (z_layer + 0.030, z_layer + 0.078):
            ax.plot([x], [y], marker='_', color='white', ms=7, mew=2.0,
                    zorder=4)

    for zc, lab in ((z_drain, r'$\phi_c^{\rm drain}$'),
                    (z_layer, r'$\phi_c^{\rm in}$')):
        ax.plot([0, 1], [zc, zc], color=CRIT, lw=1.5, zorder=5)
        ax.text(1.04, zc, '%s\n$z/H=%.2f$' % (lab, zc), fontsize=9,
                color=CRIT, va='center', ha='left')

    for i in range(1, 10):                    # the ten solved slices
        ax.plot([0, 1], [i / 10.0, i / 10.0], color='0.5', lw=0.5, ls=':',
                zorder=3)
    ax.text(-0.06, 0.5, 'ten depth slices,\none periodic cell each',
            fontsize=9, rotation=90, va='center', ha='center', color='0.35')

    ax.text(0.5, -0.025, 'cold surface, $z/H=0$', fontsize=9.5,
            ha='center', va='bottom')
    ax.text(0.5, 1.025, 'warm base, $z/H=1$', fontsize=9.5,
            ha='center', va='top')
    ax.set_xlim(-0.12, 1.55); ax.set_ylim(1.0, 0.0)
    ax.axis('off')


def cell_pocket(ax):
    ax.add_patch(Rectangle((0, 0), 1, 1, fc=ICE, ec='0.25', lw=1.2))
    rng = np.random.default_rng(1)
    for _ in range(16):
        ax.add_patch(Ellipse(rng.uniform(0.12, 0.88, 2), 0.10, 0.06,
                             angle=rng.uniform(0, 180), fc=BRINE, ec='none'))
    ax.set_title('above $\\phi_c^{\\rm in}$: isolated pockets', fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')


def cell_layered(ax):
    ax.add_patch(Rectangle((0, 0), 1, 1, fc=ICE, ec='0.25', lw=1.2))
    for x in (0.25, 0.5, 0.75):
        ax.add_patch(Rectangle((x - 0.035, 0), 0.07, 1, fc=BRINE, ec='none'))
        for y in (0.28, 0.68):               # the bridges carrying the load
            ax.add_patch(Rectangle((x - 0.035, y), 0.07, 0.09, fc=ICE,
                                   ec='0.4', lw=0.5))
    rng = np.random.default_rng(2)
    for _ in range(8):
        ax.add_patch(Ellipse(rng.uniform(0.08, 0.92, 2), 0.07, 0.045,
                             fc=BRINE, ec='none', alpha=0.8))
    ax.set_title('below $\\phi_c^{\\rm in}$: brine spans the plane', fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    z_drain, z_layer = thresholds()

    fig = plt.figure(figsize=(11.0, 5.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.15],
                          hspace=0.42, wspace=0.10)
    axc = fig.add_subplot(gs[:, 0])
    column(axc, z_drain, z_layer)
    cell_pocket(fig.add_subplot(gs[0, 1]))
    cell_layered(fig.add_subplot(gs[1, 1]))

    fig.tight_layout()
    p = os.path.join(out, 'fig_column.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)
    print('  drainage threshold at z/H = %.2f, in-plane at %.2f'
          % (z_drain, z_layer))


if __name__ == '__main__':
    main()
