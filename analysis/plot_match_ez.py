r"""Figures for the E(z) match against Gogolaze and Kujala.

(a) GOGOLAZE. Beam rigidity against the exponent, with his two reductions drawn
    as a band rather than a line: they differ by 1.81x, so the target is a
    range. The measured exponent band is shaded for comparison.

(b) KUJALA. He reports no porosity, so nothing can be matched forwards. The
    panel shows the brine profile his beams IMPLY under the closure, against
    the synthetic column. His is monotonic and reaches the skeletal range,
    ours is C-shaped.

Marchenko is deliberately absent. His profile is a measured brine profile
pushed through a correlation obtained elsewhere from three-point bending and
then fitted to a Kerr-Palmer form, so it is a construction rather than a
measurement; and its endpoint ratio disagrees with Kujala's by a factor of two
to three, so no single calibration could satisfy both. Dropping it is a
judgement about the data, and it is recorded here so the omission is not
mistaken for an oversight.

    python3 analysis/plot_match_ez.py [outdir]
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import ez_closure as ez
from match_ez import (flexural, K_TOP, K_BOT,
                      GOGO_APP, GOGO_COR, H_GOGO)
from shape_diagnosis import ours_phi

GA_CEILING = 2.0


def panel_gogolaze(ax, z):
    zc = z * H_GOGO * 100.0
    phi = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    ns = np.linspace(0.2, 5.0, 90)
    ef = np.array([flexural(ez.E_of_phi(phi, n=n, floor=0.0), z) for n in ns])
    ax.plot(ns, ef, color=fs.BLUE, lw=2.6)
    ax.axhspan(GOGO_APP, GOGO_COR, color=fs.GREEN, alpha=0.30, zorder=0)
    ax.text(3.5, 0.95, 'his two reductions\n0.785 - 1.421 GPa',
            fontsize=10, color=fs.GREEN, ha='center')
    ax.axvspan(ez.N_LO, ez.N_HI, color=fs.SKY, alpha=0.30, zorder=0)
    ax.text(ez.N_HI + 0.25, 9.0, 'measured band', fontsize=9.5,
            color=fs.BLUE, ha='left')
    for tgt, lab, c in ((GOGO_COR, 'root-corrected', fs.VERM),
                        (3.2, 'if 2.3x stiffer', fs.PURPLE)):
        try:
            nn = brentq(lambda n: flexural(ez.E_of_phi(phi, n=n, floor=0.0), z)
                        - tgt, 0.2, 12.0)
            ax.plot([nn], [tgt], 'o', color=c, ms=10, zorder=5)
            ax.annotate(r'%s: $n=%.2f$' % (lab, nn), xy=(nn, tgt),
                        xytext=(nn + 0.9, tgt * 2.4), fontsize=9.5, color=c,
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.2))
        except ValueError:
            pass
    ax.set_yscale('log'); ax.set_ylim(0.6, 14)
    ax.set_xlabel('bridge exponent $n$ in $E\\propto b^{\\,n}$')
    ax.set_ylabel(r'beam rigidity $12D/H^3$ [GPa]')
    ax.text(0.015, 0.965, '(a)', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')


def panel_kujala(ax, z):
    Et, Eb = K_TOP.mean(), K_BOT.mean()
    tgt = Et + (Eb - Et) * z
    # The closure is piecewise, so invert it branch by branch: search the
    # layered branch first and fall back to the pocket branch. A plain bracket
    # across phi_c would straddle the jump and brentq would return the
    # discontinuity rather than a root.
    for n, c, ls in ((ez.N_MID, fs.BLUE, '-'), (ez.N_LO, fs.GREEN, (0, (5, 2)))):
        req = []
        for t in tgt:
            f = lambda p: float(ez.E_of_phi(p, n=n, floor=0.0)) - t
            try:
                req.append(brentq(f, 1e-9, ez.PHI_0 - 1e-9))
            except ValueError:
                req.append(np.nan)
        ax.plot(req, z, color=c, ls=ls, lw=2.4,
                label=r'implied by Kujala, $n=%.2f$' % n)
    ax.plot(ours_phi(z), z, color=fs.VERM, lw=2.4, ls='-.',
            label='our synthetic column')
    ax.axvspan(0.25, 0.55, color=fs.ORANGE, alpha=0.16, zorder=0)
    ax.text(0.40, 0.10, 'skeletal\nrange', fontsize=10, color=fs.ORANGE,
            ha='center')

    # The ramped closure is continuous and strictly decreasing over 0..phi_0,
    # so every measured modulus has exactly one preimage and the inversion runs
    # the whole depth. The step form used earlier left a gap of moduli with no
    # solution, which cut these curves off at mid-depth.
    ax.text(0.055, 0.35, 'closure is monotone,\nso every $E$ inverts',
            fontsize=8.6, color='0.35', ha='left')
    fs.depth_axis(ax)
    ax.set_xlim(0, 0.55)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.text(0.015, 0.965, '(b)', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')
    ax.legend(loc='lower right', fontsize=9.5)


def main():
    """Two panels now. The Marchenko comparison is withdrawn with its dataset:
    his profile is a brine profile pushed through a correlation borrowed from
    three-point bending and then fitted to a functional form, so it is a
    construction rather than a measurement, and it disagrees with Kujala's
    endpoint ratio by a factor of two to three. Keeping it meant no calibration
    could satisfy the field data, because the field data do not agree."""
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    z = np.linspace(1e-3, 1.0, 400)
    fig, ax = plt.subplots(1, 2, figsize=(12.4, 5.6))
    panel_gogolaze(ax[0], z)
    panel_kujala(ax[1], z)
    fig.tight_layout()
    p = os.path.join(outdir, 'match_ez.png')
    fig.savefig(p, dpi=165)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
