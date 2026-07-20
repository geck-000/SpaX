"""Second-order (bending) size sweep on the realistic sea-ice channelled-base
microstructure, to confirm whether a genuine MCST length scale exists for the
actual morphology (prior null was for isolated spheres / generic voids).

Mirrors the validated quadratic deck rve_porous_q3.csv (L=0.24/0.32/0.40,
L_mesh=0.033, r_avg=0.04, Kappa=0.11, Mode=utx, Mode2=ss13 -> utx+ss13+ben),
but the inclusion phase is the soft brine (K/G) with a vertical channel network
at base brine fraction. Run quadratic: SPAX_MESH_ORDER=2.

Eq.19 then reads E_bending(L), E_eff (utx), G_eff (ss13) and takes the slope of
E_app vs 1/L^2 (reference-free). Expected, per prior work: slope ~ 0 -> no l.
"""
import csv
from make_ice_studies import row, E_matrix

SIZES = [0.24, 0.32, 0.40]      # L/d = 3,4,5 at d=0.08
SEEDS = 4
PHI_B = 0.13                    # warm channelled base
GAS = 0.012
CH_FRAC = 0.40

def main():
    rows = []
    E_mat = E_matrix(-4.5)
    for L in SIZES:
        for s in range(1, SEEDS + 1):
            rid = 'SI2_L%03d_s%d' % (int(L * 1000), s)
            r = row(rid, E_mat, PHI_B, GAS, 0.65, 0.65,
                    channels_frac=CH_FRAC, r_avg=0.04)
            r['L'] = '%.2f' % L
            r['r_std'] = '0.011'
            r['Mode2'] = 'Simple Shear S13'   # -> G_eff
            r['Kappa'] = '0.11'               # enable bending (-> ben deck)
            r['Bending_Plane'] = 'xz'
            r['Bending_PBC_Type'] = 'Lesicar'
            rows.append(r)
    from make_ice_studies import COLS
    with open('rve_seaice_2nd.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader(); w.writerows(rows)
    print('wrote rve_seaice_2nd.csv  (%d RVEs, %d solves: utx+ss13+ben)'
          % (len(rows), 3 * len(rows)))

if __name__ == '__main__':
    main()
