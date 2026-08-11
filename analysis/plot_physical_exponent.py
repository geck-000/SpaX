r"""Both comparisons at the exponent the geometry actually supports.

Assur's b runs 0.52-0.78 down these columns. Gibson and Ashby's square law is
derived for relative densities below about 0.3, so it is out of range at every
depth by a factor near two; at b ~ 0.5-0.7 the lamellar plane is a perforated
solid with ice in the majority, bending does not dominate, and the scaling is
nearer linear. b^1 is drawn as the defensible choice and b^0.85, what our own
cells report, beside it. b^2 is kept only as a greyed reference so the cost of
honesty is visible rather than hidden.

    python3 analysis/plot_physical_exponent.py [outdir]
"""
import os
import statistics as st
import sys
import csv
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

import layered_law as law

K_TOP = np.array([7.18, 8.16, 8.25, 8.60])
K_BOT = np.array([0.86, 1.25, 1.56, 1.42])
K_Z0 = np.array([0.37, 0.38, 0.39, 0.38])
GOGO_APP, GOGO_COR, MATRIX_FACTOR = 0.785, 1.421, 0.49
H_BEAM = 0.32
GA_VALID = 0.30

EXPS = ((1.0, r'$b^{1}$  area / stretch', fs.BLUE, '-', 2.6, True),
        (0.85, r'$b^{0.85}$  our cells', fs.GREEN, (0, (5, 2)), 2.0, True),
        (2.0, r'$b^{2}$  Gibson-Ashby (out of range)', '0.55', (0, (1, 2)),
         1.8, False))


def neutral(E, z):
    return float(np.trapz(E * z, z) / np.trapz(E, z))


def flexural(E, z):
    z0 = neutral(E, z)
    return float(12.0 * np.trapz(E * (z - z0) ** 2, z) / (z[-1] - z[0]) ** 3)


def load_col(path):
    g = defaultdict(list)
    for r in csv.DictReader(open(path, encoding='utf8', errors='replace')):
        try:
            v = float(r['E_eff'])
        except (ValueError, KeyError, TypeError):
            continue
        if v > 0:
            g[r['run_id'].split('_s')[0]].append(v / 1e9)
    ks = sorted(g, key=lambda k: int(''.join(c for c in k if c.isdigit()) or 0))
    return np.array([st.mean(g[k]) for k in ks])


def fig_kujala(outdir):
    m = load_col(os.path.join(outdir, 'results_column_ensemble.csv'))
    zc = np.linspace(0.05, 0.95, len(m))
    z = np.linspace(1e-3, 1.0, 400)
    phi = np.interp(z, [0, .29, .63, .79, .96, 1.0],
                    [0.104, 0.086, 0.128, 0.168, 0.227, 0.227])
    k = np.interp(z, zc, m)[0] / law.pocket(phi)[0]

    fig, ax = plt.subplots(1, 2, figsize=(12.8, 5.6))
    a = ax[0]
    for n, lab, c, ls, lw, inrange in EXPS:
        a.plot(law.layered(phi, n) * k, z, color=c, ls=ls, lw=lw, label=lab,
               zorder=3 if inrange else 2)
    a.plot(np.interp(z, zc, m), z, color=fs.ORANGE, lw=2.0, marker='o',
           markevery=40, ms=5, label='pocket column (current Fig. 18)')
    a.fill_betweenx([0, 1], [K_TOP.min(), K_BOT.min()],
                    [K_TOP.max(), K_BOT.max()], color=fs.VERM, alpha=0.16,
                    zorder=0, label='Kujala 1990, inferred $E(z)$')
    a.errorbar([K_TOP.mean()], [0.0], xerr=[K_TOP.std(ddof=0)], marker='^',
               ms=9, color=fs.VERM, capsize=3, zorder=5)
    a.errorbar([K_BOT.mean()], [1.0], xerr=[K_BOT.std(ddof=0)], marker='v',
               ms=9, color=fs.VERM, capsize=3, zorder=5)
    fs.depth_axis(a)
    a.set_xscale('log'); a.set_xlim(0.3, 14)
    a.set_xlabel("Effective Young's modulus [GPa]")
    a.set_title('(a) Kujala: profile at each exponent')
    a.legend(loc='lower left', fontsize=9.5, framealpha=0.95)

    b = ax[1]
    at_k = float((K_BOT / K_TOP).mean())
    theirs = [at_k, float(K_Z0.mean()),
              3 * (1 + at_k) / (at_k ** 2 + 4 * at_k + 1)]
    Ep = np.interp(z, zc, m)
    series = [('pocket', [Ep[-1] / Ep[0], neutral(Ep, z), Ep[0] / flexural(Ep, z)],
               fs.ORANGE)]
    for n, lab, c, ls, lw, inrange in EXPS:
        E = law.layered(phi, n) * k
        series.append((lab.split()[0], [E[-1] / E[0], neutral(E, z),
                                        E[0] / flexural(E, z)], c))
    series.append(('Kujala', theirs, fs.VERM))

    x = np.arange(3)
    w = 0.15
    for i, (lab, vals, c) in enumerate(series):
        b.bar(x + (i - 2) * w, vals, w, color=c, label=lab,
              alpha=0.9 if lab != r'$b^{2}$' else 0.55)
    b.set_xticks(x)
    b.set_xticklabels([r'$\alpha=E_b/E_t$', r'$z_0/H$', r'$E_t/E_{\rm flex}$'],
                      fontsize=12)
    b.set_ylabel('value')
    b.set_title('(b) the three shape metrics')
    b.legend(fontsize=9.5, ncol=2)
    b.grid(axis='x', alpha=0)

    fig.tight_layout()
    p = os.path.join(outdir, 'kujala_exponent.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)
    for lab, vals, _ in series:
        print('  %-10s alpha %.3f  z0/H %.3f  Et/Eflex %.3f'
              % (lab, *vals))


def fig_gogolaze(outdir):
    z = np.linspace(1e-3, 1.0, 600)
    zc = z * H_BEAM * 100.0
    phi = (0.29315 * zc ** 2 - 5.124 * zc + 85.977) / 1000.0

    fig, ax = plt.subplots(1, 2, figsize=(12.8, 5.4))

    a = ax[0]
    for n, lab, c, ls, lw, inrange in EXPS:
        a.plot(law.layered(phi, n), z, color=c, ls=ls, lw=lw, label=lab)
    a.plot(law.pocket(phi), z, color=fs.ORANGE, lw=2.0, label='pockets')
    a.plot(law.assur_b(phi), z, color=fs.PURPLE, lw=1.4, ls='-.',
           label=r'Assur $b$ (right axis scale)')
    a.axvline(GA_VALID, color='0.4', lw=1.2, ls=':')
    a.text(GA_VALID * 1.06, 0.10, 'Gibson-Ashby\nvalid below this $b$',
           fontsize=9.5, color='0.35')
    fs.depth_axis(a)
    a.set_xscale('log'); a.set_xlim(0.15, 14)
    a.set_xlabel("modulus [GPa]   /   $b$ [-]")
    a.set_title('(a) Gogolaze beam: $b$ never enters the valid range')
    a.legend(loc='lower left', fontsize=9.5)

    b = ax[1]
    cases = [('pockets', flexural(law.pocket(phi), z), fs.ORANGE)]
    for n, lab, c, ls, lw, inrange in EXPS:
        cases.append((lab.split()[0], flexural(law.layered(phi, n), z), c))
    cases.append((r'$\times0.49$ factor',
                  flexural(law.pocket(phi), z) * MATRIX_FACTOR, fs.PURPLE))
    xs = np.arange(len(cases))
    b.bar(xs, [c[1] for c in cases], color=[c[2] for c in cases], alpha=0.88)
    b.axhspan(GOGO_APP, GOGO_COR, color=fs.GREEN, alpha=0.32, zorder=0)
    b.text(len(cases) - 1.4, 1.03, 'measured\n0.785 - 1.421', fontsize=10,
           color=fs.GREEN, ha='center')
    for i, (_, v, _) in enumerate(cases):
        b.text(i, v * 1.08, '%.2f' % v, ha='center', fontsize=10)
    b.set_xticks(xs)
    b.set_xticklabels([c[0] for c in cases], fontsize=10.5)
    b.set_yscale('log'); b.set_ylim(0.4, 16)
    b.set_ylabel(r'flexural modulus $12D/H^3$ [GPa]')
    b.set_title('(b) whole-beam rigidity vs the fitted factor')
    b.grid(axis='x', alpha=0)

    fig.tight_layout()
    p = os.path.join(outdir, 'gogolaze_exponent.png')
    fig.savefig(p, dpi=170)
    print('wrote %s' % p)
    for lab, v, _ in cases:
        print('  %-16s %6.3f GPa   %.2fx root-corrected' % (lab, v, v / GOGO_COR))


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else '.'
    fig_kujala(outdir)
    print()
    fig_gogolaze(outdir)


if __name__ == '__main__':
    main()
