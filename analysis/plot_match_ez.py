r"""Figures for the E(z) match against Marchenko, Gogolaze and Kujala.

One panel per case, each drawn on the porosity that case supplies:

(a) MARCHENKO. His Kerr-Palmer curve against the closure evaluated on the
    porosity his own correlation implies. Drawn twice -- as computed, and
    normalised to his surface value -- because the shape agrees while the level
    does not, and the level gap is his intercept rather than our microstructure.

(b) GOGOLAZE. Beam rigidity against the exponent, with his two reductions drawn
    as a band rather than a line: they differ by 1.81x, so the target is a
    range. The Gibson-Ashby ceiling is marked, which is what the required
    exponent has to be judged against.

(c) KUJALA. He reports no porosity, so nothing can be matched forwards. The
    panel instead shows the brine profile his beams IMPLY under the closure,
    against the synthetic column we have been comparing him with. That contrast
    is the finding: his is monotonic and reaches the skeletal range, ours is
    C-shaped and stops at 0.227.

    python3 analysis/plot_match_ez.py [outdir]
"""
import os
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import layered_law as law
import ez_closure as ez
from match_ez import (marchenko_E, corr_inv, flexural, K_TOP, K_BOT,
                      GOGO_APP, GOGO_COR, H_GOGO)
from shape_diagnosis import ours_phi

GA_CEILING = 2.0


def panel_marchenko(ax, z):
    """His profile against the ADOPTED piecewise closure.

    Panels (a)-(c) all used layered_law.layered here, which applies the bridge
    factor at every depth. That is no longer the adopted form -- the closure
    branches at the percolation threshold -- so this figure and the closure
    figure were drawing two different models. All three panels now call
    ez_closure, which is what the paper states.
    """
    phi = corr_inv(marchenko_E(z))
    tgt = marchenko_E(z)

    def miss(n):
        E = ez.E_of_phi(phi, n=n, floor=0.0)
        return float(np.mean((E / E[0] - tgt / tgt[0]) ** 2))
    n = minimize_scalar(miss, bounds=(0.05, 8.0), method='bounded').x
    E = ez.E_of_phi(phi, n=ez.N_MID, floor=0.0)
    Eb = ez.E_of_phi(phi, n=n, floor=0.0)

    # where his own porosity crosses the threshold, which is where the closure
    # changes branch and where the knee sits
    zc = z[np.argmax(phi >= ez.PHI_C)] if (phi >= ez.PHI_C).any() else None

    ax.plot(tgt, z, color=fs.VERM, lw=2.6, label='Marchenko 2024 (Kerr-Palmer)')
    ax.plot(E, z, color=fs.BLUE, lw=2.2,
            label=r'closure, $n=%.2f$ (calibrated band)' % ez.N_MID)
    ax.plot(Eb, z, color=fs.PURPLE, lw=1.8, ls=(0, (4, 3)),
            label=r'best shape, $n=%.2f$ (below the band)' % n)
    if zc is not None:
        ax.axhline(zc, color='0.45', lw=1.1, ls=':')
        ax.text(0.85, zc - 0.03, r'$\phi_c$ crossed: closure changes branch',
                fontsize=8.8, color='0.35', va='bottom')
    fs.depth_axis(ax)
    ax.set_xscale('log')
    ax.set_xlim(0.3, 14)
    ax.set_xlabel("Young's modulus [GPa]")
    ax.text(0.015, 0.965, '(a)', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')
    ax.legend(loc='lower left', fontsize=8.6)
    return n


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
    ax.text(ez.N_HI + 0.25, 9.0, 'calibrated band', fontsize=9.5,
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
    ax.text(0.015, 0.965, '(b)', transform=ax.transAxes,
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
            root = np.nan
            for lo, hi in ((ez.PHI_C + 1e-6, 0.95), (1e-6, ez.PHI_C - 1e-6)):
                try:
                    if f(lo) * f(hi) < 0:
                        root = brentq(f, lo, hi)
                        break
                except ValueError:
                    pass
            req.append(root)
        ax.plot(req, z, color=c, ls=ls, lw=2.4,
                label=r'implied by Kujala, $n=%.2f$' % n)
    ax.plot(ours_phi(z), z, color=fs.VERM, lw=2.4, ls='-.',
            label='our synthetic column')
    ax.axvspan(0.25, 0.55, color=fs.ORANGE, alpha=0.16, zorder=0)
    ax.text(0.40, 0.10, 'skeletal\nrange', fontsize=10, color=fs.ORANGE,
            ha='center')

    # Where the inversion has no solution at all. The branches differ by
    # b(phi_c)^n_eff at the threshold, so a band of moduli corresponds to no
    # brine fraction. Kujala's upper column sits inside it, which is why the
    # implied-porosity curves stop rather than continue.
    n_eff = ez.N_MID * (ez.A0_REF_MM / ez.A0_MM) ** ez.SPACING_EXP
    E_hi = ez.E_ICE * (1 - 1.65 * ez.PHI_C)
    E_lo = E_hi * (1 - np.sqrt(ez.PHI_C / ez.PHI_0)) ** n_eff
    zgap = z[(tgt > E_lo) & (tgt < E_hi)]
    if len(zgap):
        ax.axhspan(zgap.min(), zgap.max(), color='0.55', alpha=0.16, zorder=0)
        ax.text(0.27, 0.5 * (zgap.min() + zgap.max()),
                'no solution:\nKujala $E$ falls in the\nbranch gap %.1f--%.1f GPa'
                % (E_lo, E_hi), fontsize=8.4, color='0.3', ha='center',
                va='center')
    fs.depth_axis(ax)
    ax.set_xlim(0, 0.55)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.text(0.015, 0.965, '(c)', transform=ax.transAxes,
            fontsize=13, fontweight='bold', va='top')
    ax.legend(loc='lower right', fontsize=9.5)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    z = np.linspace(1e-3, 1.0, 400)
    fig, ax = plt.subplots(1, 3, figsize=(17.5, 5.6))
    n = panel_marchenko(ax[0], z)
    panel_gogolaze(ax[1], z)
    panel_kujala(ax[2], z)
    fig.tight_layout()
    p = os.path.join(outdir, 'match_ez.png')
    fig.savefig(p, dpi=165)
    print('wrote %s  (Marchenko best exponent %.2f)' % (p, n))


if __name__ == '__main__':
    main()
