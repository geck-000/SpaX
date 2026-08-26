"""Read the full 6x6 Voigt stiffness per depth slice and report what the paper
needs from it.

Three things:
  1. the tensor itself at the cold top and the warm base, which section 4.2.2
     quotes;
  2. the two anisotropy ratios down the column, E_z/E_xy and G_axial/G_xy,
     which is the tensorial statement of the anisotropy result;
  3. the engineering constants per slice, which the lamination assembly of
     section 4.6 consumes (it needs E_y and G_xy, not just E_x).

Run: python coltensor.py <dir-of-elasticity_tensor_*.csv>
"""
import csv, glob, os, re, sys
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else '.'


def read_C(path):
    """The 6x6 block sits between the header line and the first blank line."""
    rows = []
    with open(path, encoding='utf8', errors='replace') as fh:
        started = False
        for line in fh:
            line = line.strip()
            if line.startswith(',11,22,33'):
                started = True
                continue
            if started:
                if not line or line.startswith('#'):
                    break
                parts = line.split(',')
                rows.append([float(x) for x in parts[1:7]])
    return np.array(rows) / 1e9        # GPa


def engineering(C):
    """Invert to the compliance and read off the engineering constants."""
    S = np.linalg.inv(C)
    E = np.array([1.0 / S[i, i] for i in range(3)])
    # file row order is 11,22,33,12,13,23 -> the shear block is G_xy, G_xz, G_yz
    G = np.array([1.0 / S[i, i] for i in range(3, 6)])   # G_xy, G_xz, G_yz
    nu = np.array([-S[0, 1] / S[0, 0], -S[0, 2] / S[0, 0], -S[1, 2] / S[1, 1]])
    return E, G, nu


files = sorted(glob.glob(os.path.join(d, 'elasticity_tensor_CTEN_z*.csv')))
if not files:
    print('no tensor files in %s' % d)
    sys.exit(1)

print('%-6s %8s %8s %8s   %8s %8s %8s   %9s %9s' %
      ('z/H', 'E_x', 'E_y', 'E_z', 'G_xy', 'G_xz', 'G_yz',
       'Ez/Exy', 'Gax/Gxy'))
rows = []
for f in files:
    z = int(re.search(r'z(\d+)', f).group(1)) / 100.0
    C = read_C(f)
    E, G, nu = engineering(C)
    Exy = 0.5 * (E[0] + E[1])
    r_E = E[2] / Exy
    r_G = 0.5 * (G[1] + G[2]) / G[0]        # (G_xz + G_yz)/2 over G_xy
    rows.append((z, E, G, r_E, r_G, C))
    print('%-6.2f %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f   %9.4f %9.4f'
          % (z, E[0], E[1], E[2], G[0], G[1], G[2], r_E, r_G))

print('\n--- the tensor at the two ends (GPa), as section 4.2.2 quotes it ---')
for label, idx in (('cold top  z/H=0.05', 0), ('warm base z/H=0.95', -1)):
    C = rows[idx][5]
    print('\n%s' % label)
    print('   C11,C22,C33 = %.3f, %.3f, %.3f' % (C[0, 0], C[1, 1], C[2, 2]))
    print('   C12,C13,C23 = %.3f, %.3f, %.3f' % (C[0, 1], C[0, 2], C[1, 2]))
    print('   C44,C55,C66 = %.3f, %.3f, %.3f' % (C[3, 3], C[4, 4], C[5, 5]))
    off = np.abs(C[:3, 3:]).max()
    print('   max |normal-shear coupling| = %.4f GPa  (%.2f%% of C11)'
          % (off, 100 * off / C[0, 0]))
    E, G, nu = engineering(C)
    print('   in-plane split E_y/E_x = %.4f' % (E[1] / E[0]))

# --- what the lamination assembly does with it -----------------------------
print('\n--- lamination assembly over the ten slices ---')
z = np.array([r[0] for r in rows])
Ex = np.array([r[1][0] for r in rows])
Ey = np.array([r[1][1] for r in rows])
nu_xy = np.array([-np.linalg.inv(r[5])[0, 1] * r[1][0] for r in rows])
h = 1.0 / len(rows)                      # equal-thickness laminae, H = 1
zc = z                                   # slice centres in normalised depth
Q = Ex / (1.0 - nu_xy ** 2)              # in-plane reduced stiffness

A = np.sum(Q * h)
# neutral plane from vanishing net axial force
z0 = np.sum(Q * h * zc) / A
D = np.sum(Q * h * ((zc - z0) ** 2 + h ** 2 / 12.0))
B = np.sum(Q * h * (zc - z0))
E_ext = A
E_flex = 12.0 * D
print('  neutral plane z0/H            = %.4f  (offset %.2f%% from mid-depth)'
      % (z0, 100 * (0.5 - z0)))
print('  E_ext = A/H                   = %.3f GPa' % E_ext)
print('  E_flex = 12D/H^3              = %.3f GPa' % E_flex)
print('  E_flex/E_ext                  = %.4f' % (E_flex / E_ext))
# coupling about the geometric mid-plane, which is what B/sqrt(AD) reports
Bm = np.sum(Q * h * (zc - 0.5))
Dm = np.sum(Q * h * ((zc - 0.5) ** 2 + h ** 2 / 12.0))
print('  B/sqrt(AD) about mid-plane    = %.4f' % (abs(Bm) / np.sqrt(A * Dm)))
print('  alpha = E_base/E_top          = %.4f' % (Ex[-1] / Ex[0]))
