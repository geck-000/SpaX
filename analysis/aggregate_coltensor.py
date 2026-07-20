#!/usr/bin/env python3
"""Aggregate the 10 per-slice full-tensor CSVs (elasticity_tensor_CTEN_z*.csv,
written by SpaX_PostProcess.extract_elasticity_tensor) into the depth
evolution of the transverse-isotropy constants + a two-panel figure.

Usage: python3 aggregate_coltensor.py <dir_with_elasticity_tensor_csvs> [out_prefix]
  -> <out_prefix>.csv  (results_coltensor.csv)  and  <out_prefix>.png

Pure numpy/matplotlib (no Abaqus). Vertical axis = z (depth); the column is
transversely isotropic about z, so the 5 TI constants are E_p (in-plane x=y),
E_z (axial), nu_p (in-plane), nu_zp (axial), G_zp (axial shear).
"""
import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VOIGT = ['11', '22', '33', '12', '13', '23']


def read_C(path):
    """Parse the 6x6 C_ij block from an elasticity_tensor_*.csv."""
    with open(path) as f:
        lines = [ln.rstrip('\n') for ln in f]
    # header row of the C block starts with ',11,22,33,12,13,23'
    hi = next(i for i, ln in enumerate(lines)
              if ln.replace(' ', '').startswith(',11,22,33'))
    C = np.zeros((6, 6))
    for r in range(6):
        parts = lines[hi + 1 + r].split(',')
        C[r, :] = [float(x) for x in parts[1:7]]
    return 0.5 * (C + C.T)   # symmetrise defensively


def eng_constants(C):
    S = np.linalg.inv(C)
    E1, E2, E3 = 1/S[0, 0], 1/S[1, 1], 1/S[2, 2]
    G12, G13, G23 = 1/S[3, 3], 1/S[4, 4], 1/S[5, 5]     # xy, xz, yz
    nu12 = -S[1, 0]*E1
    nu13 = -S[2, 0]*E1
    nu23 = -S[2, 1]*E2
    return dict(E_x=E1, E_y=E2, E_z=E3, G_xy=G12, G_xz=G13, G_yz=G23,
                nu_xy=nu12, nu_xz=nu13, nu_yz=nu23)


def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    d = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else 'results_coltensor'
    files = sorted(glob.glob(os.path.join(d, 'elasticity_tensor_CTEN_z*.csv')))
    if not files:
        print('no elasticity_tensor_CTEN_z*.csv in %s' % d); return 1

    rows = []
    for fp in files:
        rid = os.path.basename(fp).replace('elasticity_tensor_', '').replace('.csv', '')
        z = int(rid.split('_z')[1]) / 100.0
        C = read_C(fp)
        ec = eng_constants(C)
        E_p = 0.5*(ec['E_x'] + ec['E_y'])          # in-plane Young (x=y)
        G_ax = 0.5*(ec['G_xz'] + ec['G_yz'])       # axial shear (xz=yz)
        ec.update(run_id=rid, z=z,
                  E_p=E_p, E_ratio=ec['E_z']/E_p,
                  G_ax=G_ax, G_ratio=G_ax/ec['G_xy'])
        rows.append(ec)
    rows.sort(key=lambda r: r['z'])

    cols = ['run_id', 'z', 'E_x', 'E_y', 'E_z', 'E_p', 'E_ratio',
            'G_xy', 'G_xz', 'G_yz', 'G_ax', 'G_ratio',
            'nu_xy', 'nu_xz', 'nu_yz']
    out_csv = prefix + '.csv'
    with open(out_csv, 'w') as f:
        f.write(','.join(cols) + '\n')
        for r in rows:
            f.write(','.join(
                (r[c] if c == 'run_id'
                 else ('%.4f' % r[c] if c in ('z', 'E_ratio', 'G_ratio',
                                              'nu_xy', 'nu_xz', 'nu_yz')
                       else '%.4e' % r[c])) for c in cols) + '\n')
    print('wrote %s (%d slices)' % (out_csv, len(rows)))
    for r in rows:
        print('  z=%.2f  E_x=%.2f E_z=%.2f GPa  E_z/E_p=%.3f  G_ax/G_xy=%.3f'
              % (r['z'], r['E_x']/1e9, r['E_z']/1e9, r['E_ratio'], r['G_ratio']))

    # ---- figure: (a) moduli vs depth, (b) anisotropy ratios vs depth ----
    z = np.array([r['z'] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(z, [r['E_x']/1e9 for r in rows], 'o-', label=r'$E_x$')
    ax[0].plot(z, [r['E_y']/1e9 for r in rows], 's--', label=r'$E_y$')
    ax[0].plot(z, [r['E_z']/1e9 for r in rows], '^-', label=r'$E_z$')
    ax[0].set_xlabel('normalized depth $z/H$'); ax[0].set_ylabel('modulus (GPa)')
    ax[0].set_title('(a) Directional moduli vs depth')
    ax[0].legend(); ax[0].grid(alpha=0.3)

    ax[1].axhline(1.0, color='k', lw=0.8, ls=':')
    ax[1].plot(z, [r['E_ratio'] for r in rows], 'o-', label=r'$E_z/E_{xy}$ (Young)')
    ax[1].plot(z, [r['G_ratio'] for r in rows], '^-', label=r'$G_{ax}/G_{xy}$ (shear)')
    ax[1].set_xlabel('normalized depth $z/H$'); ax[1].set_ylabel('anisotropy ratio')
    ax[1].set_title('(b) Vertical anisotropy vs depth')
    ax[1].legend(); ax[1].grid(alpha=0.3)

    fig.tight_layout()
    out_png = prefix + '.png'
    fig.savefig(out_png, dpi=160)
    print('wrote %s' % out_png)
    return 0


if __name__ == '__main__':
    sys.exit(main())
