r"""Figures for the E(z) match against Gogolaze and Kujala.

(a) GOGOLAZE. Beam rigidity against the exponent, with his two reductions drawn
    as a band rather than a line: they differ by 1.81x, so the target is a
    range. The measured exponent band is shaded for comparison.

(b) KUJALA. He reports no porosity, so nothing can be matched forwards. The
    panel shows the brine profile his beams IMPLY under the closure, against
    the synthetic column. His is monotonic, 0.065 at the surface to 0.166 at
    the base, so it climbs toward the skeletal range drawn at 0.25 without
    reaching it; ours is C-shaped.

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
    over 0..phi_0 under the ramp, so every measured modulus has exactly one
    preimage; the step is not invertible, which is the one place the ramp is
    retained.
    """
    Et, Eb = K_TOP.mean(), K_BOT.mean()
    tgt = Et + (Eb - Et) * z
    out = []
    for t in tgt:
        f = lambda q: float(ez.E_of_phi(q, n=n, weight='ramp', floor=0.0)) - t
        try:
            out.append(brentq(f, 1e-9, ez.PHI_0 - 1e-9))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def panel_modulus(ax, z):
    """Both studies against one model profile."""
    lo, mid, hi = ez.E_band(ref_phi(z), floor=0.0)
    ax.fill_betweenx(z, lo, hi, color=fs.SKY, alpha=0.45, lw=0,
                     label=r'$n(b)$ band')
    ax.plot(mid, z, color=fs.BLUE, lw=1.8, label=r'closure, $n(b)$', zorder=3)

    # Kujala: moduli at two depths only, four beams each, drawn as spread.
    # clip_on=False so the surface markers are not sliced by the z=0 axis.
    for E, zz, lab in ((K_TOP, 0.015, 'Kujala, surface and base'),
                       (K_BOT, 0.985, None)):
        ax.plot([E.min(), E.max()], [zz, zz], color=fs.ORANGE, lw=1.6,
                zorder=4, clip_on=False)
        ax.plot(E, np.full_like(E, zz), 'o', color=fs.ORANGE, ms=4.2,
                mec='white', mew=0.6, label=lab, zorder=5, clip_on=False)

    # Gogolaze: a whole-beam rigidity, so it is a depth-integrated number and
    # cannot be drawn as a profile. His two reductions differ by 1.81x, which
    # is why the target is a band rather than a line. The label runs up the
    # band so it cannot collide with the legend, which sits in the empty
    # middle of the panel.
    ax.axvspan(GOGO_APP, GOGO_COR, color=fs.GREEN, alpha=0.20, lw=0, zorder=0)
    ax.text(np.sqrt(GOGO_APP * GOGO_COR), 0.40,
            'Gogolaze beam rigidity', fontsize=7.0, color='#00614A',
            ha='center', va='center', rotation=90)

    fs.depth_axis(ax)
    ax.set_ylim(1.02, -0.02)
    ax.set_xscale('log')
    ax.set_xlim(0.6, 13)
    ax.set_xticks([1, 2, 5, 10])
    ax.set_xticklabels(['1', '2', '5', '10'])
    ax.set_xlabel(r"Young's modulus $E$  (GPa)")
    fs.panel(ax, 'a')
    fs.clean(ax)
    ax.legend(loc='center', bbox_to_anchor=(0.56, 0.62))


def panel_brine(ax, z):
    """The composition each comparison rests on."""
    ax.plot(ref_phi(z), z, color=fs.BLUE, lw=1.8,
            label='our reference column')
    ax.plot(gogo_phi(z), z, color=fs.GREEN, lw=1.6, ls='-',
            label='Gogolaze, measured')
    ax.plot(kujala_implied_phi(z, None), z, color=fs.ORANGE, lw=1.6,
            ls=(0, (4.5, 1.8)), label='Kujala, implied by inversion')
    ax.axvspan(0.25, 0.55, color=fs.ORANGE, alpha=0.15, lw=0, zorder=0)
    ax.text(0.40, 0.055, 'skeletal range', fontsize=7.0, color='#8A5A00',
            ha='center', va='center')
    fs.depth_axis(ax, label=False)
    ax.set_ylim(1.02, -0.02)
    ax.set_xlim(0, 0.55)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xlabel(r'brine volume fraction $\phi$')
    fs.panel(ax, 'b')
    fs.clean(ax)
    ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.72))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    z = np.linspace(0.001, 0.999, 400)
    fig, ax = plt.subplots(1, 2, figsize=fs.size(0.46), sharey=True)
    panel_modulus(ax[0], z)
    panel_brine(ax[1], z)
    fig.tight_layout(pad=0.3, w_pad=1.0)
    p = os.path.join(outdir, 'match_ez.png')
    fig.savefig(p)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
