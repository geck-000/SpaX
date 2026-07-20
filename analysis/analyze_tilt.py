"""Study #6 analysis: does channel inclination dilute the vertical anisotropy?

Reads results_tilt{00,15,30}.csv (straight vs wavy channels leaning up to 15/30
deg, 4 seeds each) and plots the vertical/horizontal Young-modulus ratio E_z/E_x
against tilt. Produces study_tilt.png (matplotlib only).

Headline: the anisotropy is DILUTED but not removed -- E_z/E_x falls from ~1.041
(straight) to ~1.03 (tilted), i.e. the anisotropy excess drops to ~70-77% of the
vertical-channel value, saturating by ~15 deg.
"""
import csv, statistics as st
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TILTS = [(0, '00'), (15, '15'), (30, '30')]

def load(tag):
    out = []
    for r in csv.DictReader(open(f'results_tilt{tag}.csv')):
        if r.get('E_z', '') not in ('', 'ERROR', 'MISSING'):
            out.append(float(r['E_z']) / float(r['E_x']))
    return out

def main():
    data = {deg: load(tag) for deg, tag in TILTS}
    degs = [d for d, _ in TILTS]
    means = [st.mean(data[d]) for d in degs]
    sds = [st.stdev(data[d]) if len(data[d]) > 1 else 0 for d in degs]

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    # individual seeds (light) + mean±sd (dark)
    for d in degs:
        ax.scatter([d] * len(data[d]), data[d], s=34, color='#8fb8d6',
                   zorder=2, edgecolor='white', linewidth=0.5)
    ax.errorbar(degs, means, yerr=sds, fmt='o-', color='#0072B2', lw=2.2,
                ms=9, capsize=5, zorder=3, label='mean $\\pm$ s.d.')
    ax.axhline(1.0, ls='--', color='0.6', lw=1.2)
    ax.text(30, 1.001, 'isotropic', ha='right', va='bottom', color='0.5', fontsize=9)
    ax.set_xlabel('max channel tilt off vertical (deg)')
    ax.set_ylabel('vertical anisotropy  $E_z/E_x$')
    ax.set_xticks(degs)
    ax.set_title('Channel inclination dilutes the vertical anisotropy')
    # annotate the dilution
    b = means[0] - 1
    for i, d in enumerate(degs[1:], 1):
        frac = 100 * (means[i] - 1) / b
        ax.annotate(f'{frac:.0f}% of\nvertical', xy=(d, means[i]),
                    xytext=(d, means[i] - 0.006), ha='center', fontsize=8,
                    color='#22303a')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig('study_tilt.png', dpi=150)
    print('wrote study_tilt.png')
    for d in degs:
        print(f'  tilt {d:>2}deg: E_z/E_x = {st.mean(data[d]):.4f} '
              f'+/- {(st.stdev(data[d]) if len(data[d])>1 else 0):.4f}  {[round(x,3) for x in data[d]]}')

if __name__ == '__main__':
    main()
