"""Replicate (seeded) first-year C-shape column for statistical scatter envelopes.

Reviewer request (Jani #12 / Reviewer 2): the depth-profile figures show a single
packing per slice; add mean/scatter envelopes. This deck repeats the exact
first-order C-shape column physics of make_ice_studies2.study_coltensor() at every
depth slice N_SEED times, so homogenising each replicate yields the per-slice mean
and spread of the effective moduli E_x (= E_xy) and E_z. Independent packings come
for free from the per-row reseeding in SpaX_Standalone (fixed SPAX_SEED, distinct
row index -> distinct random packing), so no code change to the solver is needed.

    cd params && python3 ../studies/make_colseeds.py     # -> rve_colseeds.csv

Then, from a directory holding SpaX_Standalone.py (fixed seed = reproducible):

    SPAX_SEED=20260723 python3 SpaX_Standalone.py params/rve_colseeds.csv out_colseeds/

First-order (full_tensor=No): two load cases per RVE (Uniaxial Tension X and Z),
so 10 slices x N_SEED replicates x 2 solves.
"""
from make_ice_studies import row, write, phi_brine, E_matrix, temperature, ZS

N_SEED = 5
S_CSHAPE = [7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0]  # Cox & Weeks 1983


def study_colseeds():
    """Ten C-shape depth slices, each replicated N_SEED times (distinct run_id ->
    distinct row index -> independent packing under a fixed SPAX_SEED)."""
    rows = []
    for z, S in zip(ZS, S_CSHAPE):
        T = temperature(z, T_top=-20.0)
        phi_b = phi_brine(S, T)
        gas = 0.012 + 0.008 * (1 - z)
        sph, grw = 0.86 - 0.26 * z, 0.40 + 0.32 * z
        r_avg = 0.030 + 0.016 * z
        ch = 0.40 if phi_b > 0.05 else 0.0
        for s in range(1, N_SEED + 1):
            rid = f'CSEED_z{int(z * 100):02d}_s{s}'
            rows.append(row(rid, E_matrix(T), phi_b, gas, sph, grw,
                            channels_frac=ch, r_avg=r_avg))
    write('rve_colseeds.csv', rows)


if __name__ == '__main__':
    study_colseeds()
