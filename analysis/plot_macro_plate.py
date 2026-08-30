r"""Where the bending stiffness of the sheet actually lives.

The flexural rigidity weights each layer by its modulus AND by the square of its
lever arm, so the two compete: the base is compliant but sits at the far fibre,
the interior is stiff but sits near the neutral plane. Which wins is not obvious
and is worth drawing rather than asserting.

Drawn for both assemblies, because the answer differs between them. With the
base described as a particle composite it is stiff enough that the far fibre
dominates; described as layered it is compliant enough to drop out, and the
sheet's bending stiffness moves upward into the cold ice.

    python3 analysis/plot_macro_plate.py [outdir]
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import ez_closure as ez

COL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'results', 'results_column_ensemble.csv')

T = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])


def columns():
    phi = ez.brine_volume(T, S)
    rows = list(csv.DictReader(open(COL, encoding='utf8')))
    nu = np.array([float(r['nu_x']) for r in rows])
    Erve = np.array([float(r['E_x']) for r in rows]) / 1e9
    Ecl = np.array([float(ez.E_of_phi(p)) for p in phi])
    hyb = np.where(phi >= ez.PHI_LAYER, Ecl, Erve)
    return phi, nu, Erve, hyb


def assemble(E, nu):
    """Neutral plane and per-layer contribution to D about the mid-plane."""
    n = len(E)
    h = 1.0 / n
    z = (np.arange(n) + 0.5) * h
    Q = E / (1 - nu ** 2)
    z0 = float(np.sum(Q * h * z) / np.sum(Q * h))
    dD = Q * h * ((z - 0.5) ** 2 + h ** 2 / 12.0)   # about the geometric mid-plane
    return z, z0, dD


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    phi, nu, Erve, hyb = columns()
    z, z0p, dDp = assemble(Erve, nu)
    _, z0h, dDh = assemble(hyb, nu)

    fig, ax = plt.subplots(1, 3, figsize=fs.size(0.40), sharey=True)

    # ---- (a) the two profiles, with their neutral planes -------------------
    a = ax[0]
    a.plot(Erve, z, 'o-', color=fs.ORANGE, ms=4.4, lw=2.6,
           label='particle throughout', zorder=2)
    a.plot(hyb, z, 's-', color=fs.BLUE, ms=2.6, lw=1.2,
           label='layered below $\\phi_c$', zorder=3)
    a.axhline(0.5, color='0.55', lw=0.8, ls=':', zorder=1)
    # The two neutral planes are only 0.03 apart, so their labels are separated
    # horizontally rather than vertically: layered at the left edge, particle
    # at the right.
    for z0, c, ls in ((z0p, fs.ORANGE, '--'), (z0h, fs.BLUE, '-.')):
        a.axhline(z0, color=c, lw=1.0, ls=ls, zorder=1)
    a.text(0.68, z0h - 0.02, r'$z_0$', fontsize=8, color=fs.BLUE, va='bottom')
    a.text(13.0, z0p + 0.02, r'$z_0$', fontsize=8, color=fs.ORANGE,
           ha='right', va='top')
    a.axhspan(0.9, 1.0, color=fs.SKY, alpha=0.28, lw=0, zorder=0)
    a.text(3.2, 0.947, 'layered', fontsize=7, color=fs.BLUE, ha='center',
           va='center')
    fs.depth_axis(a)
    a.set_xscale('log'); a.set_xlim(0.6, 14)
    a.set_xticks([1, 10]); a.set_xticklabels(['1', '10'])
    a.set_xlabel(r'$E_x(z)$  (GPa)')
    fs.panel(a, 'a')
    fs.clean(a)

    # ---- (b) where the bending stiffness sits ------------------------------
    b = ax[1]
    w = 0.038
    b.barh(z - w / 2, 100 * dDp / dDp.sum(), height=w, color=fs.ORANGE,
           alpha=0.85, label='particle throughout')
    b.barh(z + w / 2, 100 * dDh / dDh.sum(), height=w, color=fs.BLUE,
           alpha=0.85, label='layered below $\\phi_c$')
    fs.depth_axis(b, label=False)
    b.set_xlabel('share of $D_{11}$  (%)')
    fs.panel(b, 'b', x=0.975, y=0.06, ha='right', va='bottom')
    fs.clean(b)
    b.grid(axis='y', alpha=0)

    # ---- (c) cumulative, read from the cold surface down -------------------
    c = ax[2]
    for dD, col, lab in ((dDp, fs.ORANGE, 'particle throughout'),
                         (dDh, fs.BLUE, 'layered below $\\phi_c$')):
        cum = 100 * np.cumsum(dD) / dD.sum()
        c.plot(cum, z, 'o-', color=col, ms=3.0, lw=1.4, label=lab)
        # how much the top half carries
        half = np.interp(0.5, z, cum)
        c.plot([half], [0.5], 'D', color=col, ms=4.5, zorder=5)
        # blue label to the right of its marker, orange to the left, so the
        # two do not stack on the mid-depth line
        if col == fs.BLUE:
            c.text(half + 5, 0.452, '%.0f%%' % half, fontsize=7.5, color=col,
                   ha='left', va='bottom')
        else:
            c.text(half - 5, 0.548, '%.0f%%' % half, fontsize=7.5, color=col,
                   ha='right', va='top')
    c.axhline(0.5, color='0.55', lw=0.8, ls=':')
    fs.depth_axis(c, label=False)
    c.set_xlim(0, 108)
    c.set_xticks([0, 25, 50, 75, 100])
    c.set_xlabel('cumulative $D_{11}$  (%)')
    fs.panel(c, 'c')
    fs.clean(c)

    # One legend for the three panels: the same two curves appear in each, so
    # repeating it three times was clutter that also cost panel width.
    h, l = a.get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.035), fontsize=8)
    fig.tight_layout(pad=0.3, w_pad=0.8, rect=(0, 0.055, 1, 1))
    p = os.path.join(out, 'fig_macro_plate.png')
    fig.savefig(p)
    print('wrote %s' % p)

    # numbers for the caption
    for lab, dD, z0 in (('particle', dDp, z0p), ('layered', dDh, z0h)):
        cum = np.cumsum(dD) / dD.sum()
        print('  %-9s z0/H %.3f | base slice %.1f%% of D | top half %.0f%%'
              % (lab, z0, 100 * dD[-1] / dD.sum(), 100 * np.interp(0.5, z, cum)))


if __name__ == '__main__':
    main()
