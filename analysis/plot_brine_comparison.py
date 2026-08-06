"""Compare brine profiles directly, instead of comparing moduli through a factor.

Marchenko (2024) does not measure a modulus profile. He measures a brine
content profile, pushes it through an empirical correlation (his formula 5,
from three-point bending tests elsewhere), and fits the result with a
Kerr & Palmer form (his Eq. 17). Inverting formula 5 on Eq. (17) therefore
recovers the brine profile his curve stands for, which can be compared with
ours without any fitted scalar at all.

(a) brine volume fraction against depth: ours as realised in the meshes,
    against the profile his curve implies.
(b) the same comparison in modulus space, showing how much of the offset is
    the correlation-versus-homogenisation difference rather than microstructure.
"""
import csv
import os
import statistics as st
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs

fs.apply()
import matplotlib.pyplot as plt

# Marchenko (2024) Eq. (17), p.11, and his empirical correlation, formula (5).
M_E0, M_ALPHA, M_N = 4.4, 0.38, 0.6
C_A, C_B = 7.23, 4.2          # E = C_A*exp(-C_B*sqrt(nu))


def marchenko(z):
    return M_E0 * (1.0 - (1.0 - M_ALPHA) * z ** M_N)


def corr(nu):
    """His formula (5): brine content -> modulus."""
    return C_A * np.exp(-C_B * np.sqrt(nu))


def corr_inv(E):
    """...and its inverse, so his curve can be read as a brine profile."""
    return (np.log(C_A / E) / C_B) ** 2


def realised(pref='MSEED', path='results_fieldseeds.csv'):
    """Per-slice realised brine fraction and effective modulus."""
    g, ph = defaultdict(list), {}
    for r in csv.DictReader(open(path)):
        v = r.get('E_eff')
        if not v or v == 'MISSING' or float(v) <= 0:
            continue
        rid = r['run_id'].rsplit('_s', 1)[0]
        if not rid.startswith(pref):
            continue
        g[rid].append(float(v) / 1e9)
        ph[rid] = float(r['phi_inclusion'])
    ks = sorted(g, key=lambda k: int(''.join(c for c in k.rsplit('_z', 1)[-1]
                                            if c.isdigit())))
    return (np.array([ph[k] for k in ks]),
            np.array([st.mean(g[k]) for k in ks]),
            np.array([st.pstdev(g[k]) for k in ks]))


def main():
    z = np.arange(0.05, 1.0, 0.1)
    Em = marchenko(z)
    nu_m = corr_inv(Em)                 # the brine profile his curve stands for
    nu_us, E_us, E_sd = realised()
    E_us_corr = corr(nu_us)             # our brine through his own correlation

    fig, ax = plt.subplots(1, 2, figsize=(12.2, 5.0))

    # ---- (a) brine against brine ----------------------------------------
    a = ax[0]
    a.plot(nu_m, z, 'k--', lw=2.0, marker='s', ms=4,
           label='implied by Marchenko 2024')
    a.plot(nu_us, z, marker='o', ms=5, color=fs.BLUE,
           label='this work, as realised')
    fs.depth_axis(a)
    a.set_xscale('log')
    a.set_xlabel('brine volume fraction')
    a.set_title('(a) the microstructure, compared without a factor')
    a.legend(loc='upper left', fontsize=10.5, framealpha=0.95)

    # ---- (b) the ratio, so the reader sees where it fails ----------------
    b = ax[1]
    ratio = nu_m / nu_us
    b.plot(ratio, z, marker='o', ms=5, color=fs.BLUE)
    b.axvline(1.0, color='0.35', lw=1.4)
    b.axvspan(0.8, 1.25, color='0.85', zorder=0)
    fs.depth_axis(b)
    b.set_xlim(0, 2.0)
    b.set_xlabel('brine implied by Marchenko / brine realised here')
    b.set_title('(b) ratio of the two profiles, by depth')
    b.annotate('basal slice: %.2f\n(we carry %.3f, he implies %.3f)'
               % (ratio[-1], nu_us[-1], nu_m[-1]),
               xy=(ratio[-1], z[-1]), xytext=(14, -2),
               textcoords='offset points', ha='left', va='center', fontsize=10)

    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig('brine_comparison.%s' % ext, dpi=200)
    print('wrote brine_comparison.{png,pdf}')

    print()
    print('%6s %12s %12s %8s' % ('z/H', 'our brine', 'his brine', 'his/ours'))
    for i, s in enumerate(z):
        print('%6.2f %12.4f %12.4f %8.2f' % (s, nu_us[i], nu_m[i], nu_m[i] / nu_us[i]))
    print()
    rr = nu_m / nu_us
    print('mean ratio his/ours over the column : %.2f' % np.mean(rr))
    print('  ... but it ranges %.2f to %.2f, so the mean is partly cancellation'
          % (rr.min(), rr.max()))
    print('  RMS departure of the ratio from unity : %.0f%%'
          % (100 * np.sqrt(np.mean((rr - 1) ** 2))))
    print('mean ratio over the nine slices above the base : %.2f' % np.mean(rr[:-1]))
    print('our brine through his correlation, vs his curve : %.1f%% RMS'
          % (100 * np.sqrt(np.mean((E_us_corr / Em - 1) ** 2))))
    print('our homogenisation vs his curve                 : %.1f%% RMS'
          % (100 * np.sqrt(np.mean((E_us / Em - 1) ** 2))))


if __name__ == '__main__':
    main()
