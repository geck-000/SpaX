r"""The sheet, the depth coordinate, and which description applies where.

Drawn rather than rendered from a mesh: the claim is about what stays the same
down the column and what changes, and no single cell can show that.

The brine is drawn where it actually sits -- inside the lamellar planes,
widening with depth and merging from separate pockets into a continuous sheet.
The plane spacing never changes, because it is set at the growth interface.

    python3 analysis/plot_column_schematic.py [outdir]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Ellipse, ConnectionPatch

import ez_closure as ez

T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])

ICE_COLD = '#eaf3fa'
ICE_WARM = '#a8c9e6'
BRINE = '#E8A33D'
BRINE_DK = '#C87F1E'
CRIT = '#C0392B'
INK = '#2b2b2b'

N_PLANES = 5


def thresholds():
    z = (np.arange(10) + 0.5) / 10.0
    phi = ez.brine_volume(T, S)
    zi = np.linspace(0, 1, 4001)
    pi = np.interp(zi, z, phi)
    return (zi[np.argmax(pi >= ez.PHI_DRAIN)],
            zi[np.argmax(pi >= ez.PHI_LAYER)], zi, pi)


def sheet(ax, zd, zl, zi, pi):
    """The sheet in vertical section."""
    cmap = LinearSegmentedColormap.from_list('ice', [ICE_COLD, ICE_WARM])
    ax.imshow(np.linspace(0, 1, 512).reshape(-1, 1), extent=[0, 1, 1, 0],
              aspect='auto', cmap=cmap, zorder=0, interpolation='bilinear')

    xs = np.linspace(0, 1, N_PLANES + 2)[1:-1]
    pitch = xs[1] - xs[0]

    # One ribbon per plane rather than a scatter of pockets. The brine occupies
    # the plane, its width follows the local brine fraction, and the eye can
    # follow a single shape from surface to base. Above the in-plane threshold
    # the ribbon is broken, because the pockets are still separate and ice
    # spans the plane; below it the ribbon is continuous and only the bridges
    # interrupt it.
    zz = np.linspace(0.012, 0.988, 600)
    half = pitch * (0.09 + 2.4 * np.interp(zz, zi, pi)) / 2.0

    for x in xs:
        m = zz >= zl
        ax.fill_betweenx(zz[m], x - half[m], x + half[m],
                         fc=BRINE, ec=BRINE_DK, lw=0.5, zorder=3)
        seg, hseg = zz[zz < zl], half[zz < zl]
        y0 = seg[0]
        while y0 < zl - 0.006:
            y1 = min(y0 + 0.020, zl - 0.006)
            k = (seg >= y0) & (seg <= y1)
            if k.sum() > 1:
                ax.fill_betweenx(seg[k], x - hseg[k], x + hseg[k],
                                 fc=BRINE, ec=BRINE_DK, lw=0.4, zorder=3)
            y0 += 0.032
        for zb in np.arange(zl + 0.030, 0.985, 0.058):
            wb = pitch * (0.09 + 2.4 * float(np.interp(zb, zi, pi)))
            ax.add_patch(Rectangle((x - wb / 2, zb - 0.009), wb, 0.018,
                                   fc=ICE_WARM, ec=BRINE_DK, lw=0.4, zorder=4))

    for zc, lab in ((zd, r'$\phi_c^{\rm drain}$,  $z/H=%.2f$' % zd),
                    (zl, r'$\phi_c^{\rm in}$,  $z/H=%.2f$' % zl)):
        ax.plot([0, 1], [zc, zc], color=CRIT, lw=1.2, zorder=6)
        ax.annotate(lab, xy=(1.0, zc), xytext=(1.09, zc),
                    fontsize=9.5, color=CRIT, va='center', ha='left',
                    arrowprops=dict(arrowstyle='-', color=CRIT, lw=1.0,
                                    shrinkA=0, shrinkB=0))

    ax.add_patch(Rectangle((0, 0), 1, 1, fc='none', ec=INK, lw=1.3, zorder=7))
    for i in range(1, 10):
        ax.plot([0, 0.038], [i / 10, i / 10], color=INK, lw=0.8, zorder=7)
        ax.plot([0.962, 1], [i / 10, i / 10], color=INK, lw=0.8, zorder=7)

    ax.annotate('', xy=(-0.10, 1.0), xytext=(-0.10, 0.0),
                arrowprops=dict(arrowstyle='->', color='0.45', lw=1.2))
    ax.text(-0.15, 0.5, 'depth  $z/H$', rotation=90, va='center', ha='center',
            fontsize=10, color='0.35')
    ax.text(0.5, -0.035, 'cold surface', ha='center', va='bottom', fontsize=10)
    ax.text(0.5, 1.035, 'warm base', ha='center', va='top', fontsize=10)
    ax.annotate('', xy=(xs[1], 0.055), xytext=(xs[2], 0.055),
                arrowprops=dict(arrowstyle='<->', color=INK, lw=1.0))
    ax.text(0.5 * (xs[1] + xs[2]), 0.085, '$a_0$', ha='center', va='top',
            fontsize=10, color=INK)
    ax.text(0.5, 1.075, 'the spacing $a_0$ is the same at every depth;\n'
            'only the brine inside a plane changes',
            ha='center', va='top', fontsize=9, color='0.3')

    ax.set_xlim(-0.20, 1.62)
    ax.set_ylim(1.20, -0.10)
    ax.axis('off')


def cell(ax, layered, title):
    ax.add_patch(Rectangle((0, 0), 1, 1, ec=INK, lw=1.2,
                           fc=ICE_WARM if layered else ICE_COLD))
    rng = np.random.default_rng(3 if layered else 7)
    if layered:
        for x in (0.28, 0.56, 0.84):
            ax.add_patch(Rectangle((x - 0.026, 0), 0.052, 1, fc=BRINE,
                                   ec=BRINE_DK, lw=0.4))
            for y in (0.30, 0.72):
                ax.add_patch(Rectangle((x - 0.026, y), 0.052, 0.072,
                                       fc=ICE_WARM, ec=BRINE_DK, lw=0.4))
        ax.annotate('ice bridge', xy=(0.306, 0.336), xytext=(0.03, 0.09),
                    fontsize=8, color=INK,
                    arrowprops=dict(arrowstyle='->', lw=0.8, color=INK))
    else:
        for x in (0.28, 0.56, 0.84):
            for y in np.arange(0.10, 1.0, 0.135):
                ax.add_patch(Ellipse((x + rng.uniform(-.012, .012),
                                      y + rng.uniform(-.02, .02)),
                                     0.085, 0.048, fc=BRINE, ec=BRINE_DK,
                                     lw=0.4))
    ax.set_title(title, fontsize=9.5, pad=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    zd, zl, zi, pi = thresholds()

    fig = plt.figure(figsize=(10.4, 5.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.30, 1.0],
                          wspace=0.04, hspace=0.34,
                          left=0.05, right=0.985, top=0.92, bottom=0.05)
    axs = fig.add_subplot(gs[:, 0])
    sheet(axs, zd, zl, zi, pi)
    a1 = fig.add_subplot(gs[0, 1])
    a2 = fig.add_subplot(gs[1, 1])
    cell(a1, False, 'above $\\phi_c^{\\rm in}$: pockets, ice spans the plane')
    cell(a2, True, 'below $\\phi_c^{\\rm in}$: brine spans it, bridges carry load')

    for ax_c, zt in ((a1, 0.42), (a2, 0.96)):
        fig.add_artist(ConnectionPatch(
            xyA=(0, 0.5), coordsA=ax_c.transAxes,
            xyB=(1.0, zt), coordsB=axs.transData,
            arrowstyle='-', color='0.6', lw=0.9, linestyle=(0, (4, 3))))

    p = os.path.join(out, 'fig_column.png')
    fig.savefig(p, dpi=180)
    print('wrote %s' % p)
    print('  drainage at z/H = %.2f, in-plane at %.2f' % (zd, zl))


if __name__ == '__main__':
    main()
