"""Re-verification of the RVE-size convergence check at the base slice.

Section 4.1.1 reports that doubling the cell edge at the base microstructure
moves the effective modulus by a few percent, and reads that as evidence the
cell is box-size converged. Two things make the original run worth repeating.

First, the results files it was drawn from predate the achieved-phase-fraction
reporting, so the comparison rested on the requested volume fractions rather
than on what the packer actually built. That distinction turned out to matter:
the achieved void fraction is now known to drift systematically with cell size
(0.050 at L/d=3 falling to 0.013 at L/d=10 against a target of 0.012 in the
size-effect sweep; the same pattern appears in results_porous_q3.csv and
results_lscale.csv). Small cells overshoot, and voids soften more than the
near-incompressible brine, so a size sweep can register a composition drift as a
size effect.

Second, the two trends point opposite ways here, which is what makes the check
informative rather than merely tidy. The void artefact makes small cells appear
SOFTER. The base sweep instead shows small cells STIFFER (4.85 GPa at L=0.50
against 4.67 at L=1.00), consistent with the channel count -- three to five
channels fit at L=0.50 against ten to eleven at L=0.80, so the percolating
network is under-represented in the small cell. If that reading is right the
achieved fractions should be flat across the sweep while the modulus falls; if
instead the achieved fractions drift, part of the reported convergence is
composition and the trend needs restating.

The microstructure is the base slice (ICE_z95) exactly as reconstructed from
results_basesweep_L065.csv, held fixed while only L varies. L_mesh is likewise
fixed at 0.033 so the inclusions are resolved identically at every size, which
is the condition that makes the comparison a size sweep rather than a mesh
sweep.

    cd params && python3 ../studies/make_basesweep.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_ice_studies import BASE, COLS

PARAMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'params')

SIZES = [0.50, 0.65, 0.80, 1.00]
N_SEED = 5

# The base slice, echoed from results_basesweep_L065.csv (which carries every
# deck field), so the re-run is the same microstructure and not a fresh guess.
MICRO = {
    'L_mesh': '0.033', 'Is_Porous': 'Composite',
    'E_matrix': '9.37e9', 'nu_matrix': '0.33',
    'VoF_sphere': '0.110', 'VoF_void_sphere': '0.020', 'VoF_incl_sphere': '0.090',
    'r_avg': '0.045', 'r_std': '0.012',
    'sphericity_avg': '0.62', 'sphericity_std': '0.10',
    'generate_channels': 'Yes', 'channel_vof_target': '0.060',
    'r_channel_avg': '0.022', 'r_channel_std': '0.005',
    'Growth_Direction': 'Z', 'Growth_Concentration': '0.70',
    'Mode': 'Uniaxial Tension X', 'Disp': '0.005',
    'Mode2': 'Uniaxial Tension Z', 'Disp2': '0.005',
    'Kappa': '0', 'full_tensor': 'No',
}


def main():
    rows = []
    for L in SIZES:
        for s in range(1, N_SEED + 1):
            r = {c: BASE.get(c, '') for c in COLS}
            r.update(MICRO)
            r['run_id'] = 'BSW_L%03d_s%d' % (round(L * 100), s)
            r['L'] = '%.2f' % L
            # Displacement is a fixed engineering strain, so it must scale with
            # the cell; a fixed 0.005 would impose a different strain at each
            # size and the sweep would no longer be a like-for-like comparison.
            r['Disp'] = r['Disp2'] = '%.4f' % (0.01 * L)
            rows.append(r)
    path = os.path.join(PARAMS, 'rve_basesweep.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print('wrote params/rve_basesweep.csv  (%d RVEs, %d solves)' % (len(rows), 2 * len(rows)))
    for L in SIZES:
        print('  L=%.2f  %d packings  elements/edge = %.0f'
              % (L, N_SEED, L / 0.033))


if __name__ == '__main__':
    main()
