#!/usr/bin/env python3
"""Study #7 -- Two-scale laminated-plate macro model of a whole sea-ice sheet.

Classical Lamination Theory (CLT) applied to the depth-resolved RVE homogenization
results: each of the 10 column slices is one graded layer with its own in-plane
(plane-stress) stiffness Q(z). Integrating Q through the thickness gives the sheet
extensional [A], bending--extension coupling [B], and bending [D] stiffness
matrices. Because the warm base is soft, the sheet is graded -> the geometric
mid-plane is NOT the neutral plane, so B != 0 (measurable bend--stretch coupling),
and the effective flexural modulus is much less than the extensional modulus.

Pure post-processing of results_column.csv (no Abaqus). Effective moduli in the
CSV use the high-frequency (~9.4 GPa) matrix; pass --beam to rescale to the
vibrating-beam-effective matrix (x0.49) for field comparison.

Usage: python3 macro_plate.py results_column.csv [H_metres] [--beam]
"""
import sys, csv
import numpy as np

BEAM_FACTOR = 0.49


def q_plane_stress(Ex, Ey, nu_xy, Gxy):
    """Plane-stress reduced stiffness (orthotropic layer), Voigt 1,2,6."""
    nu_yx = nu_xy * Ey / Ex
    d = 1.0 - nu_xy * nu_yx
    Q = np.zeros((3, 3))
    Q[0, 0] = Ex / d
    Q[1, 1] = Ey / d
    Q[0, 1] = Q[1, 0] = nu_xy * Ey / d
    Q[2, 2] = Gxy
    return Q


def main():
    if len(sys.argv) < 2:
        print(__doc__); return 1
    csv_path = sys.argv[1]
    H = 1.0
    beam = '--beam' in sys.argv
    for a in sys.argv[2:]:
        if a != '--beam':
            H = float(a)
    scale = BEAM_FACTOR if beam else 1.0

    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    rows.sort(key=lambda r: r['run_id'])
    n = len(rows)
    t = H / n                                   # equal-thickness layers
    # layer interfaces z_k, top = -H/2, base = +H/2 (slice z05 near top)
    zk = np.linspace(-H / 2.0, H / 2.0, n + 1)

    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    prof = []
    for k, r in enumerate(rows):
        Ex = float(r['E_x']) * scale
        Ey = float(r['E_y']) * scale
        Gxy = float(r['G_xy']) * scale
        nu_xy = float(r.get('nu_x', r.get('nu_eff', 0.33)))
        Q = q_plane_stress(Ex, Ey, nu_xy, Gxy)
        z1, z2 = zk[k], zk[k + 1]
        A += Q * (z2 - z1)
        B += Q * 0.5 * (z2**2 - z1**2)
        D += Q * (1.0 / 3.0) * (z2**3 - z1**3)
        zc = 0.5 * (z1 + z2)
        prof.append((r['run_id'], (k + 0.5) / n, Ex, zc, Q[0, 0]))

    # Effective engineering moduli of the sheet
    # extensional: E_ext = A11*(1 - nu^2)/H  (approx, from A^-1)
    Ainv = np.linalg.inv(A)
    E_ext = 1.0 / (Ainv[0, 0] * H)
    # flexural (bending) effective modulus: E_flex = 12/H^3 / Dinv11-equivalent
    Dinv = np.linalg.inv(D)
    E_flex = 12.0 / (Dinv[0, 0] * H**3)
    # neutral-plane offset from geometric mid-plane (bend-stretch coupling)
    z_na = B[0, 0] / A[0, 0]
    coupling = abs(B[0, 0]) / np.sqrt(A[0, 0] * D[0, 0])   # dimensionless B/sqrt(AD)

    print("=" * 66)
    print("LAMINATED-PLATE MACRO MODEL  (%d layers, H=%.3f m, matrix=%s)"
          % (n, H, 'beam-eff x0.49' if beam else 'high-freq'))
    print("=" * 66)
    print("  A11 (extensional) = %.4e N/m" % A[0, 0])
    print("  D11 (bending)     = %.4e N*m" % D[0, 0])
    print("  B11 (coupling)    = %.4e N   (0 only if ungraded)" % B[0, 0])
    print("  neutral plane offset from mid-plane: %.4f * H (= %.1f%% toward cold top)"
          % (z_na / H, 100 * z_na / H))
    print("  dimensionless bend-stretch coupling B/sqrt(AD) = %.4f" % coupling)
    print("  effective EXTENSIONAL modulus E_ext = %.3f GPa" % (E_ext / 1e9))
    print("  effective FLEXURAL   modulus E_flex = %.3f GPa" % (E_flex / 1e9))
    print("  flexural/extensional = %.3f (soft base knocks bending down more)"
          % (E_flex / E_ext))

    out = 'results_macro_plate.csv'
    with open(out, 'w') as f:
        f.write('# Laminated-plate macro model from %s (H=%.3f m, %s matrix)\n'
                % (csv_path, H, 'beam' if beam else 'highfreq'))
        f.write('quantity,value,unit\n')
        f.write('A11,%.6e,N/m\n' % A[0, 0])
        f.write('B11,%.6e,N\n' % B[0, 0])
        f.write('D11,%.6e,N*m\n' % D[0, 0])
        f.write('neutral_plane_offset_over_H,%.6f,-\n' % (z_na / H))
        f.write('coupling_B_over_sqrtAD,%.6f,-\n' % coupling)
        f.write('E_extensional,%.6e,Pa\n' % E_ext)
        f.write('E_flexural,%.6e,Pa\n' % E_flex)
        f.write('E_flex_over_E_ext,%.6f,-\n' % (E_flex / E_ext))
        f.write('\n# per-layer profile\nrun_id,z_over_H,E_x_Pa,z_centroid_m,Q11_Pa\n')
        for rid, zoh, Ex, zc, Q11 in prof:
            f.write('%s,%.3f,%.6e,%.6f,%.6e\n' % (rid, zoh, Ex, zc, Q11))
    print("\n  wrote %s" % out)

    # ---- figure: (a) in-plane modulus vs depth; (b) bending-weight Q*z^2 vs depth
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    zoh = np.array([p[1] for p in prof])
    Exs = np.array([p[2] for p in prof]) / 1e9
    zc = np.array([p[3] for p in prof])
    Q11 = np.array([p[4] for p in prof])
    bend_contrib = Q11 * zc**2 * t                  # contribution of each layer to D11
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(zoh, Exs, 'o-')
    ax[0].axvline((0.5 - z_na / H), color='r', ls='--', lw=1,
                  label='neutral plane')
    ax[0].set_xlabel('normalized depth $z/H$'); ax[0].set_ylabel('$E_x$ (GPa)')
    ax[0].set_title('(a) In-plane modulus vs depth')
    ax[0].legend(); ax[0].grid(alpha=0.3)
    ax[1].bar(zoh, bend_contrib / bend_contrib.sum() * 100, width=0.8 / len(zoh) * 1.0)
    ax[1].set_xlabel('normalized depth $z/H$')
    ax[1].set_ylabel('% of bending stiffness $D_{11}$')
    ax[1].set_title('(b) Where the bending stiffness lives')
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig('study_macro_plate.png', dpi=160)
    print("  wrote study_macro_plate.png")
    return 0


if __name__ == '__main__':
    sys.exit(main())
