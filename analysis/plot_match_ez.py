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
from match_ez import (marchenko_E, corr_inv, flexural, K_TOP, K_BOT,
                      GOGO_APP, GOGO_COR, H_GOGO)
from shape_diagnosis import ours_phi

GA_CEILING = 2.0


def panel_marchenko(ax, z):
    phi = corr_inv(marchenko_E(z))
    tgt = marchenko_E(z)

    def miss(n):
        E = law.layered(phi, n)
        return float(np.mean((E / E[0] - tgt / tgt[0]) ** 2))
    n = minimize_scalar(miss, bounds=(0.2, 8.0), method='bounded').x
    E = law.layered(phi, n)

    ax.plot(tgt, z, color=fs.VERM, lw=2.6, label='Marchenko 2024 (Kerr-Palmer)')
    ax.plot(E, z, color=fs.BLUE, lw=2.2,
            label=r'closure, $b^{%.2f}$ (as computed)' % n)
    ax.plot(E * tgt[0] / E[0], z, color=fs.BLUE, lw=1.8, ls=(0, (4, 3)),
            label='same, normalised to his surface')
    ax.annotate('', xy=(tgt[0], 0.06), xytext=(E[0], 0.06),
                arrowprops=dict(arrowstyle='<->', color='0.35', lw=1.3))
    ax.text(4.9, 0.13, r'$\times%.2f$: his intercept' % (E[0] / tgt[0]),
            fontsize=10, color='0.35')
    fs.depth_axis(ax)
    ax.set_xlim(0.8, 11)
    ax.set_xlabel("Young's modulus [GPa]")
    ax.set_title('(a) Marchenko: shape matches, level is his intercept')
    ax.legend(loc='lower left', fontsize=9.5)
    return n


def panel_gogolaze(ax, z):
    zc = z * H_GOGO * 100.0
    phi = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    ns = np.linspace(0.4, 5.0, 90)
    ef = np.array([flexural(law.layered(phi, n), z) for n in ns])
    ax.plot(ns, ef, color=fs.BLUE, lw=2.6)
    ax.axhspan(GOGO_APP, GOGO_COR, color=fs.GREEN, alpha=0.30, zorder=0)
    ax.text(3.5, 0.95, 'his two reductions\n0.785 - 1.421 GPa',
            fontsize=10, color=fs.GREEN, ha='center')
    ax.axvspan(0.4, GA_CEILING, color=fs.SKY, alpha=0.18, zorder=0)
    ax.text(1.15, 8.5, 'cellular-solid\nrange', fontsize=10, color=fs.BLUE,
            ha='center')
    for tgt, lab, c in ((GOGO_COR, 'root-corrected', fs.VERM),
                        (3.2, 'if 2.3x stiffer', fs.PURPLE)):
        try:
            nn = brentq(lambda n: flexural(law.layered(phi, n), z) - tgt,
                        0.3, 12.0)
            ax.plot([nn], [tgt], 'o', color=c, ms=10, zorder=5)
            ax.annotate(r'%s: $b^{%.2f}$' % (lab, nn), xy=(nn, tgt),
                        xytext=(nn - 1.5, tgt * 2.1), fontsize=10, color=c,
                        arrowprops=dict(arrowstyle='->', color=c, lw=1.2))
        except ValueError:
            pass
    ax.set_yscale('log'); ax.set_ylim(0.6, 14)
    ax.set_xlabel('bridge exponent $n$ in $E\\propto b^{\\,n}$')
    ax.set_ylabel(r'beam rigidity $12D/H^3$ [GPa]')
    ax.set_title('(b) Gogolaze: required exponent vs what theory allows')


def panel_kujala(ax, z):
    Et, Eb = K_TOP.mean(), K_BOT.mean()
    tgt = Et + (Eb - Et) * z
    for n, c, ls in ((1.0, fs.BLUE, '-'), (2.0, fs.GREEN, (0, (5, 2)))):
        req = []
        for t in tgt:
            try:
                req.append(brentq(lambda p: law.layered(p, n) - t, 1e-6, 0.95))
            except ValueError:
                req.append(np.nan)
        ax.plot(req, z, color=c, ls=ls, lw=2.4,
                label=r'implied by Kujala, $b^{%.0f}$' % n)
    ax.plot(ours_phi(z), z, color=fs.VERM, lw=2.4, ls='-.',
            label='our synthetic column')
    ax.axvspan(0.25, 0.55, color=fs.ORANGE, alpha=0.16, zorder=0)
    ax.text(0.40, 0.10, 'skeletal\nrange', fontsize=10, color=fs.ORANGE,
            ha='center')
    fs.depth_axis(ax)
    ax.set_xlim(0, 0.55)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.set_title('(c) Kujala: the porosity his beams imply')
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
