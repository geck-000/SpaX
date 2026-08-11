r"""Figures for the layered basal knockdown and the drained/undrained bracket.

Four panels, in the order the argument runs:

(a) the plate-spacing law, which is the mechanism -- transverse modulus against
    layer pitch, drained and undrained, with the physical spacing of sea ice
    marked. The two limits diverge as the spacing falls, which is what makes
    drainage state decisive for a layered morphology and irrelevant for pockets.
(b) the column profile, pocket against layered, with Kujala's measured surface
    and base moduli. The neutral-axis position is annotated because it is the
    quantity fixed by shape alone.
(c) Gogolaze's cantilever as flexural modulus, where the bracket has to contain
    a whole-beam value rather than a profile.
(d) where the layered cells land against the empirical E(v_b) correlations,
    which is the argument in one frame: the pocket law is far shallower than
    every field correlation, and the drained/undrained bracket spans them. The
    undrained points saturate at the ice modulus at low porosity, which is both
    physically right -- sealed brine in a thin layer resists at its bulk
    modulus -- and the ceiling of the multiplicative combination used here.

Run from results/:  python3 ../analysis/plot_layered_bracket.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

from gogolaze_layered import (brine, pocket, layered, blend, flexural,
                              H_BEAM, E_MEAS_APPARENT, E_MEAS_CORRECTED,
                              MATRIX_FACTOR)

# layer-count sweep, phi = 0.15, L = 0.5, cluster campaign bracket_nlayers
N = np.array([1, 2, 3, 4])
A0 = 0.5 / N
E_DRN = np.array([0.942, 0.704, 0.490, 0.357])
E_UND = np.array([2.640, 4.782, 5.052, 5.442])
SD_DRN = np.array([0.044, 0.014, 0.010, 0.015])
SD_UND = np.array([0.412, 0.384, 0.039, 0.007])

# Kujala et al. (1990), four strain-gauged beams
K_TOP, K_BOT = np.array([7.18, 8.16, 8.25, 8.60]), np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = (0.37, 0.39)

CELL_MM = (2.0, 5.0)      # model unit cell edge in mm, from the paper
PLATE_MM = (0.5, 1.0)     # sea-ice plate spacing


def panel_spacing(ax):
    ax.errorbar(A0, E_UND, yerr=SD_UND, marker='o', color=fs.VERM,
                label='undrained (brine sealed)', capsize=3)
    ax.errorbar(A0, E_DRN, yerr=SD_DRN, marker='s', color=fs.BLUE,
                label='drained (pressure relaxed)', capsize=3)
    ax.fill_between(A0, E_DRN, E_UND, color=fs.SKY, alpha=0.15, zorder=0)
    # physical plate spacing, converted into model units via the cell edge
    lo = PLATE_MM[0] / CELL_MM[1] * 0.5
    hi = PLATE_MM[1] / CELL_MM[0] * 0.5
    ax.axvspan(lo, hi, color=fs.GREEN, alpha=0.15, zorder=0)
    ax.axhspan(K_BOT.min(), K_BOT.max(), color=fs.ORANGE, alpha=0.22, zorder=0)
    ax.text(0.135, 1.9, 'sea-ice\nplate spacing', fontsize=11,
            color=fs.GREEN, ha='center')
    ax.text(0.42, 1.05, 'Kujala base', fontsize=11, color=fs.ORANGE)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'layer pitch $a_0$  (model units)')
    ax.set_ylabel(r'$E_x$  (GPa)')
    ax.set_title('(a) the mechanism: plate spacing sets the modulus')
    ax.legend(loc='lower left', fontsize=11)


def panel_profile(ax):
    z = np.linspace(0, 1, 400)
    phi = np.interp(z, [0, .29, .63, .79, .96, 1.0],
                    [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])
    Ep = pocket(phi)
    w = np.clip((z - 0.75) / 0.25, 0.0, 1.0)
    El = blend(w, phi, True)
    for E, c, ls, lw, lab in ((Ep, fs.VERM, '--', 3.4, 'pockets throughout'),
                              (El, fs.BLUE, '-', 2.2, 'layered base, drained')):
        z0 = np.trapz(E * z, z) / np.trapz(E, z)
        ax.plot(E, z, color=c, ls=ls, lw=lw,
                label=r'%s  ($z_0/H=%.2f$)' % (lab, z0))
        ax.axhline(z0, color=c, ls=':', lw=1.4)
    ax.errorbar(K_TOP.mean(), 0.0, xerr=K_TOP.std(), marker='D',
                color=fs.BLACK, capsize=4, label='Kujala surface / base')
    ax.errorbar(K_BOT.mean(), 1.0, xerr=K_BOT.std(), marker='D',
                color=fs.BLACK, capsize=4)
    ax.axhspan(K_Z0[0], K_Z0[1], color=fs.GREEN, alpha=0.20, zorder=0)
    ax.text(6.4, 0.345, r'measured $z_0/H$', fontsize=11, color=fs.GREEN)
    fs.depth_axis(ax)
    ax.set_xscale('log'); ax.set_xlim(0.1, 12)
    ax.set_xlabel(r"Young's modulus  (GPa)")
    ax.set_title('(b) the column, against Kujala (1990)')
    ax.legend(loc='lower left', fontsize=10.5)


def panel_gogolaze(ax):
    z = np.linspace(0, 1, 600)
    phi = brine(z * H_BEAM * 100.0)
    one = np.ones_like(z)
    cases = [
        ('pockets\n(current)', flexural(pocket(phi), z)[0], fs.VERM),
        (r'pockets $\times0.49$' '\n(fudge factor)',
         flexural(pocket(phi) * MATRIX_FACTOR, z)[0], fs.PURPLE),
        ('layers\nundrained', flexural(blend(one, phi, False), z)[0], fs.ORANGE),
        ('layers\ndrained', flexural(blend(one, phi, True), z)[0], fs.BLUE),
    ]
    xs = np.arange(len(cases))
    ax.bar(xs, [c[1] for c in cases], color=[c[2] for c in cases], alpha=0.85)
    ax.axhspan(E_MEAS_APPARENT, E_MEAS_CORRECTED, color=fs.GREEN, alpha=0.30,
               zorder=0)
    ax.text(2.4, 1.05, 'Gogolaze measured\n0.785 - 1.421 GPa', fontsize=11,
            color=fs.GREEN, ha='center')
    for x, (_, v, _) in zip(xs, cases):
        ax.text(x, v * 1.10, '%.2f' % v, ha='center', fontsize=11)
    ax.set_xticks(xs); ax.set_xticklabels([c[0] for c in cases], fontsize=11)
    ax.set_yscale('log'); ax.set_ylim(0.15, 20)
    ax.set_ylabel(r'flexural modulus $12D/H^3$  (GPa)')
    ax.set_title('(c) Gogolaze cantilever: bracket contains the beam')
    ax.grid(axis='x', alpha=0)


def panel_laws(ax):
    """Where the layered cells land against the empirical E(v_b) correlations.

    The pocket law is far shallower than every field correlation, which is the
    gap the whole exercise is about. Plotting the layered cells on the same axes
    shows whether the mechanism reaches them, without any fitted scale.
    """
    v = np.linspace(0.005, 0.26, 300)
    ax.plot(v, 9.5 * (1 - np.sqrt(v)) ** 4, color=fs.BLACK,
            label='Weeks & Assur 1967')
    ax.plot(v, 7.23 * np.exp(-4.2 * np.sqrt(v)), color=fs.PURPLE, ls='--',
            label='Marchenko 2024')
    ax.plot(v, 3.1031 * np.exp(-3.385 * np.sqrt(v)), color=fs.GREEN, ls=':',
            label='Karulina 2019')
    ax.plot(v, 9.37 * (1 - 1.65 * v), color=fs.VERM, lw=2.6,
            label='SpaX pockets (this work)')

    phi_pts = np.array([0.10, 0.15, 0.227])
    ax.plot(phi_pts, layered(phi_pts, True), 's', color=fs.BLUE, ms=9,
            label='SpaX layered, drained')
    ax.plot(phi_pts, layered(phi_pts, False), 'o', color=fs.ORANGE, ms=9,
            mfc='none', mew=2, label='SpaX layered, undrained')
    for p in phi_pts:
        ax.plot([p, p], [layered(np.array([p]), True)[0],
                         layered(np.array([p]), False)[0]],
                color=fs.SKY, lw=6, alpha=0.35, zorder=0)

    ax.set_yscale('log'); ax.set_ylim(0.1, 12)
    ax.set_xlabel(r'brine volume fraction $v_b$')
    ax.set_ylabel(r'$E$  (GPa)')
    ax.set_title('(d) layered cells reach the empirical correlations')
    ax.legend(fontsize=10, loc='lower left', ncol=2)


def main():
    fig, ax = plt.subplots(2, 2, figsize=(13.5, 11))
    panel_spacing(ax[0, 0])
    panel_profile(ax[0, 1])
    panel_gogolaze(ax[1, 0])
    panel_laws(ax[1, 1])
    fig.tight_layout()
    out = sys.argv[1] if len(sys.argv) > 1 else 'layered_bracket.png'
    fig.savefig(out, dpi=170)
    print('wrote %s' % out)


if __name__ == '__main__':
    main()
