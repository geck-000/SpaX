r"""Figures for the usable E(z) closure.

(a) E(phi) with the exponent band, against the empirical correlations and our
    own pocket cells. The pocket law is the cold-end limit the closure returns
    as b tends to one, so the two curves meeting there is a check, not a fit.
(b) E(z) for a stated column, with the band, against Kujala's measured surface
    and base moduli. This is the practical output.
(c) the two datasets the exponent was calibrated on, each on its own porosity,
    showing what the band does and does not cover.
(d) sensitivity: how far each ingredient moves the answer, with the measured
    ones separated from the assumed and calibrated ones. This is the panel that
    says how much to trust the rest.

    python3 analysis/plot_ez_closure.py [outdir]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import ez_closure as ez

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
GOGO = (0.785, 1.421)


def panel_ephi(ax):
    phi = np.linspace(0.002, 0.24, 400)
    lo, mid, hi = ez.E_band(phi, floor=0.0)
    ax.fill_between(phi, lo, hi, color=fs.SKY, alpha=0.35,
                    label=r'closure, $n=%.2f$-$%.2f$' % (ez.N_LO, ez.N_HI))
    ax.plot(phi, mid, color=fs.BLUE, lw=2.6, label=r'closure, $n=%.2f$' % ez.N_MID)
    ax.plot(phi, ez.E_ICE * (1 - 1.65 * phi), color=fs.ORANGE, lw=2.0, ls='--',
            label='pocket cells (measured)')
    ax.plot(phi, 9.5 * (1 - np.sqrt(phi)) ** 4, color=fs.BLACK, lw=1.6,
            label='Weeks & Assur 1967')
    ax.plot(phi, 7.23 * np.exp(-4.2 * np.sqrt(phi)), color=fs.PURPLE, lw=1.6,
            ls=':', label='Marchenko 2024')
    ax.axvline(ez.PHI_0, color='0.45', lw=1.2, ls='-.')
    ax.text(ez.PHI_0 * 1.02, 6.2, r'$\phi_0$: plane empties', fontsize=9.5,
            color='0.35', rotation=90, va='top')
    ax.set_yscale('log'); ax.set_ylim(0.05, 12); ax.set_xlim(0, 0.24)
    ax.set_xlabel(r'brine volume fraction $\phi$')
    ax.set_ylabel(r'$E$  [GPa]')
    ax.set_title('(a) the closure against measured laws')
    ax.legend(fontsize=9, loc='lower left')


def panel_column(ax):
    z = np.linspace(0, 1, 300)
    for n, c, lw, ls, lab in ((ez.N_LO, fs.SKY, 1.6, '--', None),
                              (ez.N_MID, fs.BLUE, 2.6, '-',
                               r'$n=%.2f$' % ez.N_MID),
                              (ez.N_HI, fs.SKY, 1.6, '--', None)):
        ax.plot(ez.E_column(z, -20.0, -1.8, 6.0, n=n), z, color=c, lw=lw,
                ls=ls, label=lab)
    lo = ez.E_column(z, -20.0, -1.8, 6.0, n=ez.N_HI)
    hi = ez.E_column(z, -20.0, -1.8, 6.0, n=ez.N_LO)
    ax.fill_betweenx(z, lo, hi, color=fs.SKY, alpha=0.30, label='exponent band')
    ax.errorbar([K_TOP.mean()], [0.0], xerr=[K_TOP.std(ddof=0)], marker='^',
                ms=10, color=fs.VERM, capsize=4, label='Kujala surface / base')
    ax.errorbar([K_BOT.mean()], [1.0], xerr=[K_BOT.std(ddof=0)], marker='v',
                ms=10, color=fs.VERM, capsize=4)
    fs.depth_axis(ax)
    ax.set_xscale('log'); ax.set_xlim(0.15, 14)
    ax.set_xlabel(r"Young's modulus  [GPa]")
    ax.set_title(r'(b) $E(z)$: 1 m column, $-20\,^\circ$C, $S=6$ ppt')
    ax.legend(fontsize=9.5, loc='lower left')


def panel_cases(ax):
    z = np.linspace(1e-3, 1.0, 400)
    zc = z * 32.0
    phi_g = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0
    fl = [ez.flexural(ez.E_of_phi(phi_g, n=n), z)
          for n in (ez.N_HI, ez.N_MID, ez.N_LO)]
    ax.bar([0], [fl[2] - fl[0]], 0.5, bottom=[fl[0]], color=fs.BLUE, alpha=0.7,
           label='closure band')
    ax.plot([0], [fl[1]], 'o', color=fs.BLUE, ms=11, zorder=5)
    ax.bar([1], [GOGO[1] - GOGO[0]], 0.5, bottom=[GOGO[0]], color=fs.GREEN,
           alpha=0.7, label="Gogolaze's two reductions")
    tgt = 4.4 * (1 - 0.62 * z ** 0.6)
    phi_m = (np.log(7.23 / tgt) / 4.2) ** 2
    al = [float(ez.E_of_phi(phi_m, n=n, floor=0.0)[-1] /
                ez.E_of_phi(phi_m, n=n, floor=0.0)[0])
          for n in (ez.N_HI, ez.N_MID, ez.N_LO)]
    ax2 = ax.twinx()
    ax2.bar([2.6], [al[2] - al[0]], 0.5, bottom=[al[0]], color=fs.BLUE,
            alpha=0.7)
    ax2.plot([2.6], [al[1]], 'o', color=fs.BLUE, ms=11, zorder=5)
    ax2.bar([3.6], [0.02], 0.5, bottom=[0.384 - 0.01], color=fs.VERM, alpha=0.8)
    ax2.set_ylabel(r'grading $\alpha=E_b/E_t$'); ax2.set_ylim(0, 0.55)
    ax2.grid(False)
    ax.set_xticks([0, 1, 2.6, 3.6])
    ax.set_xticklabels(['closure', 'measured', 'closure', 'measured'],
                       fontsize=10)
    ax.set_xlim(-0.6, 4.3)
    ax.set_ylabel(r'beam rigidity  [GPa]'); ax.set_ylim(0, 2.4)
    ax.set_title('(c) Gogolaze level (left) and Marchenko grading (right)')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(axis='x', alpha=0)


def panel_sensitivity(ax):
    phi_ref = 0.12
    base = float(ez.E_of_phi(phi_ref))
    rows = [
        ('exponent $n$\n0.49-0.59  CALIBRATED',
         float(ez.E_of_phi(phi_ref, n=ez.N_LO)),
         float(ez.E_of_phi(phi_ref, n=ez.N_HI)), fs.ORANGE),
        (r'$\phi_0$  0.15-0.36' '\nASSUMED',
         float(ez.E_of_phi(phi_ref, phi_0=0.36)),
         float(ez.E_of_phi(phi_ref, phi_0=0.15)), fs.VERM),
        (r'$a_0$  0.20-0.50 mm' '\nMEASURED',
         float(ez.E_of_phi(phi_ref, a0_mm=0.50)),
         float(ez.E_of_phi(phi_ref, a0_mm=0.20)), fs.GREEN),
    ]
    for i, (lab, lo, hi, c) in enumerate(rows):
        ax.barh(i, hi - lo, left=lo, height=0.45, color=c, alpha=0.85)
        ax.text(max(hi, lo) * 1.04, i, r'$\times$%.1f' % (max(hi, lo) / min(hi, lo)),
                va='center', fontsize=10.5)
    ax.axvline(base, color=fs.BLACK, lw=1.6, ls='--')
    ax.text(base * 1.04, -0.62, 'nominal %.2f GPa' % base, fontsize=10,
            ha='left', va='bottom')
    ax.set_ylim(-0.8, len(rows) - 0.4)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
    ax.set_xscale('log'); ax.set_xlim(0.2, 6)
    ax.set_xlabel(r'$E$ at $\phi=%.2f$  [GPa]' % phi_ref)
    ax.set_title('(d) what each ingredient is worth')
    ax.grid(axis='y', alpha=0)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    fig, ax = plt.subplots(2, 2, figsize=(13.6, 11))
    panel_ephi(ax[0, 0])
    panel_column(ax[0, 1])
    panel_cases(ax[1, 0])
    panel_sensitivity(ax[1, 1])
    fig.tight_layout()
    p = os.path.join(outdir, 'ez_closure.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)


if __name__ == '__main__':
    main()
