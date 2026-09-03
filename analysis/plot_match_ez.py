r"""Figures for the E(z) match against Gogolaze and Kujala.

(a) GOGOLAZE. Beam rigidity against the exponent, with his two reductions drawn
    as a band rather than a line: they differ by 1.81x, so the target is a
    range. The measured exponent band is shaded for comparison.

(c) LANGLEBEN AND POUNDER. Their acoustic fits are the only field comparison
    made at HIGH FREQUENCY, so they measure the instantaneous modulus the
    closure computes and need no level correction. Both of their fits stop at
    phi = 0.10, which on the reference column is z/H = 0.889, so the panel
    draws them solid where they are supported and dashed below that.

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


# --- Langleben & Pounder (1963), from the primary source -------------------
# Their Eqs. 8 and 9, converted from 10^10 dyn/cm^-2 (= 1 GPa) and from brine
# volume as a PERCENTAGE to the fraction phi used here:
#   Eq. 8  cold annual ice   E = (10.00 - 0.351 v) -> 10.00 - 35.1 phi   GPa
#   Eq. 9  warm polar ice    E = ( 8.90 - 0.163 v) ->  8.90 - 16.3 phi   GPa
# Both fitted over v <= 10 %, the range Figs. 1 and 2 of that chapter span.
LP_COLD = (10.00, 35.1)
LP_WARM = (8.90, 16.3)
LP_PHI_MAX = 0.10


def lp_E(phi, coeff):
    a, b = coeff
    return a - b * np.asarray(phi, float)


def panel_langleben(ax, z):
    """The closure against the one high-frequency dataset, in E(phi(z)).

    This panel exists because the comparison it draws needs NO level
    correction. A resonance measurement is fast enough that the delayed elastic
    strain has not developed, so what Langleben and Pounder report is the
    instantaneous modulus -- the same quantity the closure returns. Kujala's
    and Gogolaze's whole-beam numbers are low-frequency and do need the
    correction, which is why they sit in panel (a) and this is separate.

    Drawn on the same reference column and the same logarithmic modulus axis as
    panel (a), so the two can be read against each other directly.
    """
    phi = ref_phi(z)
    lo, mid, hi = ez.E_band(phi, floor=0.0)
    ax.fill_betweenx(z, lo, hi, color=fs.SKY, alpha=0.45, lw=0)
    ax.plot(mid, z, color=fs.BLUE, lw=1.8, label=r'closure, $n(b)$', zorder=3)

    # supported where phi <= 0.10, extrapolated below; the split is the point.
    ok = phi <= LP_PHI_MAX
    z_split = z[ok].max() if ok.any() else np.nan
    for coeff, col, lab in ((LP_COLD, fs.VERM, 'L&P cold annual ice'),
                            (LP_WARM, fs.PURPLE, 'L&P warm polar ice')):
        E = lp_E(phi, coeff)
        ax.plot(np.where(ok, E, np.nan), z, color=col, lw=1.6, label=lab,
                zorder=4)
        ax.plot(np.where(~ok, E, np.nan), z, color=col, lw=1.2,
                ls=(0, (2.2, 1.8)), zorder=4)

    if np.isfinite(z_split):
        ax.axhline(z_split, color=fs.GREY, lw=0.8, ls=':', zorder=1)
        ax.text(0.68, z_split - 0.018,
                r'$\phi=0.10$: limit of both fits', fontsize=6.8,
                color='#4A4A4A', ha='left', va='bottom')

    fs.depth_axis(ax, label=False)
    ax.set_ylim(1.02, -0.02)
    ax.set_xscale('log')
    ax.set_xlim(0.6, 13)
    ax.set_xticks([1, 2, 5, 10])
    ax.set_xticklabels(['1', '2', '5', '10'])
    ax.set_xlabel(r"Young's modulus $E$  (GPa)")
    fs.panel(ax, 'c')
    fs.clean(ax)
    ax.legend(loc='center left', bbox_to_anchor=(0.02, 0.60), fontsize=6.8)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    z = np.linspace(0.001, 0.999, 400)
    fig, ax = plt.subplots(1, 3, figsize=fs.size(0.32), sharey=True)
    panel_modulus(ax[0], z)
    panel_brine(ax[1], z)
    panel_langleben(ax[2], z)
    fig.tight_layout(pad=0.3, w_pad=0.9)
    p = os.path.join(outdir, 'match_ez.png')
    fig.savefig(p)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
