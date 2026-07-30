"""Nonlocal length-scale study: extend the three-size bending sweep.

Why
---
The production sweep (rve_seaice_2nd.csv -> results_si2nd.csv) bends the
channelled-base RVE at L/d = 3, 4, 5 with four packings each, and the apparent
bending modulus *rises* with cell size (0.911, 0.937, 1.000 of the classical
plate modulus). Read through the modified couple-stress form

    E_app/E_inf = 1 + 12 l^2 / L^2                       (type III.A)

this gives a negative slope, hence l^2 < 0 and no couple-stress length scale.
That is the result reported in the manuscript, and it is correct as far as it
goes -- but it tests only one of the two nonclassical families. Srinivasa &
Reddy (Appl. Mech. Rev. 69(3) 030802, 2017) classify gradient/couple-stress
models (III.A, stiffening as the body shrinks) and integral nonlocal models of
Eringen type (III.B, softening) as formally dual, related by exchanging the
roles of stress and strain. Softening is precisely what we measure, so the same
data read through

    E_inf/E_app = 1 + (e0a)^2 / L^2                      (type III.B)

returns a *real* nonlocal length. Fitting the three published points gives
e0a = 1.05 d with an intercept of 0.964 where it should be 1.000 -- suggestive,
but three points cannot pin an asymptote.

What this deck adds
-------------------
1. Two larger cells (L/d = 6, 8). The intercept is the large-cell asymptote, so
   it is constrained almost entirely by the widest cells; the existing sweep
   has none beyond L/d=5. L/d=10 is included at reduced seed count as a
   stretch: the L=0.96 cells of the earlier bending studies failed to solve
   (results_bending.csv and results_lscale.csv both carry zeros there), so
   L=0.80 is treated as the practical ceiling and anything above it as a bonus.
2. Six packings per size instead of four, since the quantity being fitted is a
   ~8% trend across the sweep while the packing scatter is a few percent.
3. A homogeneous baseline at *matched* cell sizes (rve_eringen_homog.csv,
   VoF=0). This is the control the existing calibration lacks: a homogeneous
   cube has no microstructure and therefore no length scale of any kind, so
   whatever size dependence it shows is the cube-versus-plate kinematics and
   the discretisation, not the material. The existing four-point calibration
   (results_homog.csv) sits a constant +12% above the classical plate modulus
   but varies only -1.0% over a 2x size range, against the +8.1% microstructural
   trend and in the opposite direction -- so the baseline does not explain the
   effect and subtracting it should slightly strengthen it. Matched sizes let
   that subtraction be done point-by-point rather than by interpolation.

Each RVE solves utx (first-order reference E_eff) and, through Kappa>0, the
bending deck -- so the ratio is formed per cell rather than against a single
global plate modulus, which removes the matrix-modulus and Poisson bookkeeping
from the fit.

Run the bending decks with SPAX_MESH_ORDER=2: linear tetrahedra lock in
bending, and the whole measurement is a few-percent modulus difference.
"""
import csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import COLS

PARAMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'params')

D = 0.08                     # mean inclusion diameter, 2 * r_avg
SIZES = [(3, 6), (4, 6), (5, 6), (6, 6), (8, 6), (10, 3)]   # (L/d, n_seeds)

# the channelled-base microstructure of rve_seaice_2nd.csv, held fixed
MICRO = {
    'L_mesh': '0.033', 'Is_Porous': 'Composite',
    'E_matrix': '9.3922e+09', 'nu_matrix': '0.33',
    'VoF_sphere': '0.09', 'r_avg': '0.040', 'r_std': '0.011',
    'Mode': 'Uniaxial Tension X', 'Disp': '0.005',
    'Mode2': 'Simple Shear S13', 'Disp2': '0.005',
    'VoF_void_sphere': '0.012', 'VoF_incl_sphere': '0.078',
    'E_sphere_inclusion': '2.2e9', 'nu_sphere_inclusion': '0.48',
    'sphericity_avg': '0.65', 'sphericity_std': '0.10',
    'min_distance': '0.002', 'max_iterations': '200000',
    'nlgeom_flag': 'OFF', 'PBC_Method': 'Gmsh',
    'Kappa': '0.11', 'Bending_Plane': 'xz', 'Bending_PBC_Type': 'Lesicar',
    'generate_channels': 'Yes', 'channel_vof_target': '0.052',
    'r_channel_avg': '0.022', 'r_channel_std': '0.005',
    'Growth_Direction': 'Z', 'Growth_Concentration': '0.65',
    'Inclusion_Type': 'Liquid', 'K_inclusion': '2200000000.0',
    'G_inclusion': '440029.33528897085', 'full_tensor': 'No',
}


def write(name, rows):
    path = os.path.join(PARAMS, name)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote params/%s  (%d RVEs, %d solves)' % (name, len(rows), 2 * len(rows)))


def study_channelled():
    rows = []
    for ratio, nseed in SIZES:
        L = round(ratio * D, 3)
        for s in range(1, nseed + 1):
            r = dict(MICRO)
            r['run_id'] = 'ERG_L%03d_s%d' % (round(L * 1000), s)
            r['L'] = '%.2f' % L
            rows.append(r)
    write('rve_eringen.csv', rows)


def study_homogeneous():
    """Same cells with the microstructure removed: the geometric baseline."""
    rows = []
    for ratio, _ in SIZES:
        L = round(ratio * D, 3)
        r = dict(MICRO)
        r.update({
            'run_id': 'ERGH_L%03d' % round(L * 1000), 'L': '%.2f' % L,
            'VoF_sphere': '0.0', 'VoF_void_sphere': '0.0',
            'VoF_incl_sphere': '0.0',
            'generate_channels': 'No', 'channel_vof_target': '0',
        })
        rows.append(r)
    write('rve_eringen_homog.csv', rows)


if __name__ == '__main__':
    study_channelled()
    study_homogeneous()
    print('\nL/d ->  L      (d = %.2f)' % D)
    for ratio, n in SIZES:
        print('  %2d   %.2f   %d packings' % (ratio, ratio * D, n))
