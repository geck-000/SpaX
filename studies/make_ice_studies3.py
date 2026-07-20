"""Third battery of sea-ice column studies (2026-06-29 pm).
Reuses helpers from make_ice_studies.py.

  9  orient -> rve_orient.csv  (crystal/pocket orientation + columnar vs granular)
  10 gas    -> rve_gas.csv     (air-void fraction sweep at fixed brine)
  11 seeds  -> rve_seeds.csv   (statistical replicates: 3 configs x 5 random packings)

Note on replicates: SpaX_Standalone seeds each row by its ROW INDEX
(seed = SPAX_SEED + idx*2654435761), so N identical-parameter rows with distinct
run_ids automatically pack differently -- that is the replicate ensemble.
"""
from make_ice_studies import row, write, E_matrix

# ---------------------------------------------------------------------------
# Study 9: CRYSTAL / POCKET ORIENTATION.
#  Part A - isolate pocket elongation axis: fixed phi_b=0.06, NO channels,
#           elongated pockets (sphericity 0.60), vary Growth_Direction (the
#           ellipsoid major-axis) X vs Z and concentration (alignment strength)
#           vs Random. Tests whether crystal-pocket orientation alone makes
#           E_z != E_x.
#  Part B - columnar (vertical channels, strong Z texture) vs granular (random,
#           round, no channels) at phi_b=0.08: the real columnar/granular contrast.
def _set_dir(r, d):
    r['Growth_Direction'] = d
    return r

def study_orient():
    gas = 0.012
    E_mat = E_matrix(-8.0)
    rows = []
    # Part A: pocket-orientation, no channels, elongated
    phiA, sphA = 0.06, 0.60
    rows.append(_set_dir(row('ORI_rand', E_mat, phiA, gas, sphA, 0.00), 'Random'))
    rows.append(_set_dir(row('ORI_Z50', E_mat, phiA, gas, sphA, 0.50), 'Z'))
    rows.append(_set_dir(row('ORI_Z90', E_mat, phiA, gas, sphA, 0.90), 'Z'))
    rows.append(_set_dir(row('ORI_X50', E_mat, phiA, gas, sphA, 0.50), 'X'))
    rows.append(_set_dir(row('ORI_X90', E_mat, phiA, gas, sphA, 0.90), 'X'))
    # Part B: columnar vs granular
    rows.append(_set_dir(row('ORI_columnar', E_mat, 0.08, gas, 0.65, 0.90,
                             channels_frac=0.40), 'Z'))
    rows.append(_set_dir(row('ORI_granular', E_mat, 0.08, gas, 0.85, 0.00), 'Random'))
    write('rve_orient.csv', rows)

# ---------------------------------------------------------------------------
# Study 10: GAS / POROSITY sweep. Fixed small brine (phi_b=0.02 meshed pockets),
# sweep the true-void air fraction 0..0.10. Isolates the gas phase's stiffness
# knockdown (gas = hard void; brine = soft solid -- both near-void for E, but gas
# carries no pressure at all). Relevant to drained/summer & gassy ice.
def study_gas():
    E_mat = E_matrix(-10.0)
    rows = []
    for g in (0.00, 0.02, 0.04, 0.06, 0.08, 0.10):
        rows.append(row(f'GAS_v{int(g*100):02d}', E_mat, 0.02, g, 0.80, 0.40))
    write('rve_gas.csv', rows)

# ---------------------------------------------------------------------------
# Study 11: STATISTICAL REPLICATES. 3 representative configs x 5 packings each
# -> scatter (std / CoV) on E and E_z/E_x, separating real signal from
# realization noise. top = cold isolated, mid = near-percolation, base = warm
# channelled.
def study_seeds():
    rows = []
    configs = [
        ('top',  E_matrix(-17.0), 0.02, 0.015, 0.84, 0.45, 0.0),
        ('mid',  E_matrix(-10.0), 0.05, 0.012, 0.75, 0.55, 0.0),
        ('base', E_matrix(-4.5),  0.08, 0.013, 0.65, 0.65, 0.40),
    ]
    for tag, E_mat, phi_b, gas, sph, grw, ch in configs:
        for s in range(1, 6):
            rows.append(row(f'SEED_{tag}_s{s}', E_mat, phi_b, gas, sph, grw,
                            channels_frac=ch))
    write('rve_seeds.csv', rows)

if __name__ == '__main__':
    study_orient()
    study_gas()
    study_seeds()
