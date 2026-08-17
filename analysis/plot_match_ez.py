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


T_REF = np.array([-19.1, -17.3, -15.4, -13.6, -11.8, -10.0, -8.2, -6.3, -4.5, -2.7])
S_REF = np.array([7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0])


def ref_phi(z):
    """The reference column of Section 4.5, interpolated onto z.

    The same profile Figure 10(b) and the sheet assembly use. Keeping one
    column across the three figures matters: the closure is monotone in phi, so
    a different column is a different sheet, not a different view of one.
    """
    zc = (np.arange(len(T_REF)) + 0.5) / len(T_REF)
    return np.interp(z, zc, ez.brine_volume(T_REF, S_REF))


def gogo_phi(z):
    """Gogolaze's measured brine profile, his eq. (14), on normalised depth."""
    zc = z * H_GOGO * 100.0
    return (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0


def kujala_implied_phi(z, n):
    """The brine profile Kujala's beams require under the closure.

    He reports no composition, so nothing can be matched forwards. Inverting
    his moduli is the only way to use him, and it asks a fair question: is the
    phi(z) his beams imply a physically plausible one? The closure is monotone
    over 0..phi_0, so every measured modulus has exactly one preimage.
    """
    Et, Eb = K_TOP.mean(), K_BOT.mean()
    tgt = Et + (Eb - Et) * z
    out = []
    for t in tgt:
        f = lambda q: float(ez.E_of_phi(q, n=n, floor=0.0)) - t
        try:
            out.append(brentq(f, 1e-9, ez.PHI_0 - 1e-9))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def panel_modulus(ax, z):
    """Both studies against one model profile."""
    lo, mid, hi = (ez.E_of_phi(ref_phi(z), n=n, floor=0.0)
                   for n in (ez.N_HI, ez.N_MID, ez.N_LO))
    ax.fill_betweenx(z, lo, hi, color=fs.SKY, alpha=0.35, label='exponent band')
    ax.plot(mid, z, color=fs.BLUE, lw=2.6, label=r'closure, $n=%.2f$' % ez.N_MID)

    # Kujala: moduli at two depths only, four beams each, drawn as spread.
    for E, zz, lab in ((K_TOP, 0.015, 'Kujala, surface and base'),
                       (K_BOT, 0.985, None)):
        ax.plot(E, np.full_like(E, zz), 'o', color=fs.ORANGE, ms=7,
                mec='white', mew=0.8, label=lab, zorder=5)
        ax.plot([E.min(), E.max()], [zz, zz], color=fs.ORANGE, lw=2.2, zorder=4)

    # Gogolaze: a whole-beam rigidity, so it is a depth-integrated number and
    # cannot be drawn as a profile. His two reductions differ by 1.81x, which
    # is why the target is a band rather than a line.
    ax.axvspan(GOGO_APP, GOGO_COR, color=fs.GREEN, alpha=0.18, zorder=0)
    ax.text(0.5 * (GOGO_APP + GOGO_COR), 0.52,
            "Gogolaze\nbeam rigidity\n(two reductions)", fontsize=8.6,
            color=fs.GREEN, ha='center', va='center')

    fs.depth_axis(ax)
    ax.set_xscale('log')
    ax.set_xlabel("Young's modulus   [GPa]")
    ax.text(0.015, 0.965, '(a)', transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='top')
    ax.legend(loc='center left', fontsize=9)


def panel_brine(ax, z):
    """The composition each comparison rests on."""
    ax.plot(ref_phi(z), z, color=fs.BLUE, lw=2.6,
            label='our reference column')
    ax.plot(gogo_phi(z), z, color=fs.GREEN, lw=2.4, ls='-',
            label='Gogolaze, measured')
    ax.plot(kujala_implied_phi(z, ez.N_MID), z, color=fs.ORANGE, lw=2.4,
            ls=(0, (5, 2)), label='Kujala, implied by inversion')
    ax.axvspan(0.25, 0.55, color=fs.ORANGE, alpha=0.14, zorder=0)
    ax.text(0.40, 0.12, 'skeletal\nrange', fontsize=9.5, color=fs.ORANGE,
            ha='center')
    fs.depth_axis(ax)
    ax.set_xlim(0, 0.55)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.text(0.015, 0.965, '(b)', transform=ax.transAxes, fontsize=13,
            fontweight='bold', va='top')
    ax.legend(loc='lower right', fontsize=9)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    z = np.linspace(0.001, 0.999, 400)
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.4))
    panel_modulus(ax[0], z)
    panel_brine(ax[1], z)
    fig.tight_layout()
    p = os.path.join(outdir, 'match_ez.png')
    fig.savefig(p, dpi=165)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
