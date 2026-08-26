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

    fig, ax = plt.subplots(1, 3, figsize=(15.2, 5.4))

    # ---- (a) the two profiles, with their neutral planes -------------------
    a = ax[0]
    a.plot(Erve, z, 'o-', color=fs.ORANGE, ms=6, lw=2.0,
           label='particle throughout')
    a.plot(hyb, z, 's-', color=fs.BLUE, ms=6, lw=2.4,
           label='layered below $\\phi_c$')
    for z0, c, ls in ((z0p, fs.ORANGE, '--'), (z0h, fs.BLUE, '-.')):
        a.axhline(z0, color=c, lw=1.3, ls=ls)
    a.axhline(0.5, color='0.5', lw=1.0, ls=':')
    a.text(11.5, 0.494, 'mid-depth', fontsize=8, color='0.45', ha='right')
    a.text(0.70, z0h - 0.014, '$z_0$ layered', fontsize=8.5, color=fs.BLUE)
    a.text(11.5, z0p + 0.040, '$z_0$ particle', fontsize=8.5, color=fs.ORANGE,
           ha='right')
    a.axhspan(0.9, 1.0, color=fs.SKY, alpha=0.25, zorder=0)
    a.text(4.6, 0.945, 'layered', fontsize=9, color=fs.BLUE, ha='center')
    fs.depth_axis(a)
    a.set_xscale('log'); a.set_xlim(0.6, 14)
    a.set_xlabel(r'$E_x(z)$   [GPa]')
    a.text(0.015, 0.965, '(a)', transform=a.transAxes, fontsize=13,
           fontweight='bold', va='top')
    a.legend(loc='lower left', fontsize=9)

    # ---- (b) where the bending stiffness sits ------------------------------
    b = ax[1]
    w = 0.038
    b.barh(z - w / 2, 100 * dDp / dDp.sum(), height=w, color=fs.ORANGE,
           alpha=0.85, label='particle throughout')
    b.barh(z + w / 2, 100 * dDh / dDh.sum(), height=w, color=fs.BLUE,
           alpha=0.85, label='layered below $\\phi_c$')
    fs.depth_axis(b)
    b.set_xlabel('share of $D_{11}$   [%]')
    b.text(0.015, 0.965, '(b)', transform=b.transAxes, fontsize=13,
           fontweight='bold', va='top')
    b.legend(loc='center right', fontsize=9)
    b.grid(axis='y', alpha=0)

    # ---- (c) cumulative, read from the cold surface down -------------------
    c = ax[2]
    for dD, col, lab in ((dDp, fs.ORANGE, 'particle throughout'),
                         (dDh, fs.BLUE, 'layered below $\\phi_c$')):
        cum = 100 * np.cumsum(dD) / dD.sum()
        c.plot(cum, z, 'o-', color=col, ms=5, lw=2.2, label=lab)
        # how much the top half carries
        half = np.interp(0.5, z, cum)
        c.plot([half], [0.5], 'D', color=col, ms=9, zorder=5)
        c.text(half + 3, 0.455 if col == fs.BLUE else 0.545,
               '%.0f%%' % half, fontsize=10, color=col)
    c.axhline(0.5, color='0.5', lw=1.0, ls=':')
    fs.depth_axis(c)
    c.set_xlim(0, 105)
    c.set_xlabel('cumulative share of $D_{11}$   [%]')
    c.text(0.015, 0.965, '(c)', transform=c.transAxes, fontsize=13,
           fontweight='bold', va='top')
    c.legend(loc='lower right', fontsize=9)

    fig.tight_layout()
    p = os.path.join(out, 'fig_macro_plate.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)

    # numbers for the caption
    for lab, dD, z0 in (('particle', dDp, z0p), ('layered', dDh, z0h)):
        cum = np.cumsum(dD) / dD.sum()
        print('  %-9s z0/H %.3f | base slice %.1f%% of D | top half %.0f%%'
              % (lab, z0, 100 * dD[-1] / dD.sum(), 100 * np.interp(0.5, z, cum)))


if __name__ == '__main__':
    main()
