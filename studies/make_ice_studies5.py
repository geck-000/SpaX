"""Fifth battery of sea-ice studies (2026-07-05).

  #3 sizechan   -> rve_sizechan.csv    (channelled morphology, RVE-size convergence
                                        of E_z/E_x: 5 box sizes x 4 seeds)
  #4 salfamily  -> rve_salfamily.csv   (salinity-profile family: measured-Arctic +
                                        linear S(z), mapped to E(z) and anisotropy)

Both are first-order (utx + utz -> E_x, E_z), linear tets. Reuses the physics
helpers from make_ice_studies.py. Seeds are distinct because SpaX_Standalone
derives the packing RNG from SPAX_SEED + row-index.
"""
from make_ice_studies import row, write, phi_brine, E_matrix, temperature, ZS

# ---------------------------------------------------------------------------
# Study 3: RVE-SIZE CONVERGENCE for the CHANNELLED morphology (first-order E).
# The generic size study used spheres; repeat for the percolated vertical channel
# network to confirm E_z/E_x is box-size-converged (not a finite-cell artifact).
# Fixed well-percolated base morphology (phi_b=0.08, channels carry 40% of brine);
# sweep box size L at fixed inclusion/channel sizes -> larger box samples more of
# the network. 4 seeds per size for scatter.
def study_sizechan():
    sizes = (0.30, 0.40, 0.50, 0.65, 0.80)
    seeds = (1, 2, 3, 4)
    E_mat = E_matrix(-6.0)                 # warm, well-percolated base
    phi_b, gas, sph, grw = 0.08, 0.012, 0.65, 0.65
    rows = []
    for L in sizes:
        for s in seeds:
            rid = f'SZCH_L{int(L*100):02d}_s{s}'
            r = row(rid, E_mat, phi_b, gas, sph, grw, channels_frac=0.40, r_avg=0.035)
            r['L'] = f'{L:.2f}'            # override box size (row() bakes the default)
            rows.append(r)
    write('rve_sizechan.csv', rows)

# ---------------------------------------------------------------------------
# Study 4: SALINITY-PROFILE FAMILY. Existing set = C-shape / monotonic / steeper.
# Add (a) a measured-Arctic-station FY profile and (b) a linear profile, and map
# how S(z) shape sets E(z) and the base anisotropy. High-frequency matrix (trend
# convention, comparable to the column/percolation studies).
def study_salfamily():
    # (a) representative measured Arctic FY station: high surface, interior minimum,
    #     strong basal rise (distinct from the idealized C-shape).
    S_arctic = [6.5, 5.0, 4.2, 3.8, 3.6, 3.7, 4.1, 4.8, 6.2, 9.0]
    # (b) linear increase with depth 3 -> 9 ppt.
    S_linear = [3.0 + 6.0 * z for z in ZS]
    rows = []
    for tag, prof in (('ARC', S_arctic), ('LIN', S_linear)):
        for z, S in zip(ZS, prof):
            T = temperature(z, T_top=-20.0)
            phi_b = min(phi_brine(S, T), 0.16)      # cap for mesh safety
            gas = 0.012
            sph, grw = 0.86 - 0.26*z, 0.40 + 0.32*z
            r_avg = 0.030 + 0.016*z
            ch = 0.40 if phi_b > 0.05 else 0.0
            rows.append(row(f'SAL_{tag}_z{int(z*100):02d}', E_matrix(T), phi_b, gas,
                            sph, grw, channels_frac=ch, r_avg=r_avg))
    write('rve_salfamily.csv', rows)

if __name__ == '__main__':
    study_sizechan()
    study_salfamily()
