"""Study #5 analysis: temperature-dependent brine bulk modulus K(T).

Reads the paired first-order columns results_brineKconst.csv (fixed K=2.2 GPa) and
results_brineKtemp.csv (physically-varying K(T)), which share byte-identical
morphology per slice, and quantifies the pure K(T) sensitivity of the effective
moduli E_x(z), E_z(z) and the anisotropy E_z/E_x(z). Produces study_brineK.png
(pandas + matplotlib only, no Abaqus).

Headline: the brine is near-incompressible, so a physical K(T) span (~1.24x down
the column) moves the effective sea-ice moduli by << 1% -- the column response is
robust to the brine's thermal-compaction (concentration) state.
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ZS = [0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]
# K(T) applied per slice (GPa), from make_ice_studies6.K_brine on the C-shape column
K_GPA = [2.777,2.752,2.718,2.676,2.625,2.567,2.499,2.423,2.339,2.247]

def load(path):
    out = {}
    for r in csv.DictReader(open(path)):
        z = r['run_id'].split('_z')[1]
        out[z] = {k: float(r[k]) for k in ('E_x','E_z','E_anisotropy')}
    return out

def main():
    kc = load('results_brineKconst.csv')
    kt = load('results_brineKtemp.csv')
    zk = [f'{int(z*100):02d}' for z in ZS]
    Exc = [kc[z]['E_x']/1e9 for z in zk]; Ext = [kt[z]['E_x']/1e9 for z in zk]
    Ezc = [kc[z]['E_z']/1e9 for z in zk]; Ezt = [kt[z]['E_z']/1e9 for z in zk]
    dEx = [100*(kt[z]['E_x']-kc[z]['E_x'])/kc[z]['E_x'] for z in zk]
    dEz = [100*(kt[z]['E_z']-kc[z]['E_z'])/kc[z]['E_z'] for z in zk]
    anc = [kc[z]['E_anisotropy'] for z in zk]; ant = [kt[z]['E_anisotropy'] for z in zk]
    mx = max(max(abs(d) for d in dEx), max(abs(d) for d in dEz))

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.4))
    # (a) applied K(T) profile
    ax[0].plot(K_GPA, ZS, 'o-', color='#1f77b4')
    ax[0].axvline(2.2, ls='--', color='0.6', lw=1, label='fixed K = 2.2 GPa')
    ax[0].set_xlabel('brine bulk modulus $K$ (GPa)'); ax[0].set_ylabel('depth $z/H$')
    ax[0].invert_yaxis(); ax[0].set_title('(a) applied $K(T)$'); ax[0].legend(fontsize=8)
    # (b) E(z): const vs temp (nearly overlapping)
    ax[1].plot(Exc, ZS, 'o-',  color='#1f77b4', label='$E_x$, fixed $K$')
    ax[1].plot(Ext, ZS, 'x--', color='#7fc7ff', label='$E_x$, $K(T)$')
    ax[1].plot(Ezc, ZS, 's-',  color='#d62728', label='$E_z$, fixed $K$')
    ax[1].plot(Ezt, ZS, '+--', color='#ff9896', label='$E_z$, $K(T)$')
    ax[1].set_xlabel('effective Young modulus (GPa)'); ax[1].set_ylabel('depth $z/H$')
    ax[1].invert_yaxis(); ax[1].set_title('(b) $E(z)$: curves overlie'); ax[1].legend(fontsize=7)
    # (c) percent shift from K(T)
    ax[2].plot(dEx, ZS, 'o-', color='#1f77b4', label='$\\Delta E_x$')
    ax[2].plot(dEz, ZS, 's-', color='#d62728', label='$\\Delta E_z$')
    ax[2].axvline(0, color='0.6', lw=1)
    ax[2].set_xlabel('shift from $K(T)$ (%)'); ax[2].set_ylabel('depth $z/H$')
    ax[2].set_xlim(-0.2, 0.2); ax[2].invert_yaxis()
    ax[2].set_title(f'(c) uniform $+$stiffening, $\\leq {mx:.2f}\\%$'); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig('study_brineK.png', dpi=140)
    print('wrote study_brineK.png')

    print(f"\n{'z/H':>5} {'K(GPa)':>7} {'dEx%':>7} {'dEz%':>7} "
          f"{'aniso_const':>11} {'aniso_temp':>10}")
    for i, z in enumerate(zk):
        print(f"{ZS[i]:>5.2f} {K_GPA[i]:>7.3f} {dEx[i]:>7.3f} {dEz[i]:>7.3f} "
              f"{anc[i]:>11.4f} {ant[i]:>10.4f}")
    print(f"\nMax |dE| across the whole column: {mx:.3f} %  "
          f"(K span 2.247->2.777 GPa, x1.236)")

if __name__ == '__main__':
    main()
