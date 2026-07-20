"""Second battery of sea-ice column studies (2026-06-29 pm).
Reuses the physics helpers from make_ice_studies.py.

  5 channel  -> rve_channel.csv     (fixed phi_b, sweep channel radius x channel fraction)
  6 fymy     -> rve_fymy.csv        (first-year C-shape column vs desalinated multi-year)
  7 basetensor -> rve_basetensor.csv (channelled base slices, full 6x6 tensor)
  8 mono2    -> rve_mono2.csv       (steeper monotonic salinity -> closer Marchenko shape)
"""
from make_ice_studies import (row, write, phi_brine, E_matrix, temperature,
                              ZS, BEAM_FACTOR)

# ---------------------------------------------------------------------------
# Study 5: CHANNEL GEOMETRY -- fixed phi_b=0.08 (well percolated), vary the
# vertical channel network: its radius and the share of brine it carries.
# Extends the percolation finding (connectivity drives E_z/E_x).
def study_channel():
    phi_b, gas = 0.08, 0.012
    E_mat = E_matrix(-6.0)
    sph, grw = 0.65, 0.65
    rows = []
    for rch in (0.015, 0.020, 0.025, 0.030):
        for frac in (0.2, 0.4, 0.6):
            rid = f'CHAN_r{int(rch*1000):03d}f{int(frac*100):02d}'
            r = row(rid, E_mat, phi_b, gas, sph, grw, channels_frac=frac)
            r['r_channel_avg'] = f'{rch:.3f}'
            r['r_channel_std'] = f'{rch*0.25:.4f}'
            rows.append(r)
    write('rve_channel.csv', rows)

# ---------------------------------------------------------------------------
# Study 6: FY vs MULTI-YEAR ice. Same thermal profile (T_top=-20), but
#   FY: first-year C-shape salinity, low gas (brine-dominated).
#   MY: desalinated (~1.5 ppt, drained over summer) with high gas (drained
#       brine channels become air) -> low brine, gassy. Is_Porous kept Composite;
#       gas is the VoF_void_sphere true-void phase.
def study_fymy():
    S_cshape = [7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0]
    rows = []
    # first-year
    for z, S in zip(ZS, S_cshape):
        T = temperature(z, T_top=-20.0)
        phi_b = phi_brine(S, T)
        gas = 0.012 + 0.008 * (1 - z)
        sph, grw = 0.86 - 0.26*z, 0.40 + 0.32*z
        r_avg = 0.030 + 0.016*z
        ch = 0.40 if phi_b > 0.05 else 0.0
        rows.append(row(f'FYMY_fy_z{int(z*100):02d}', E_matrix(T), phi_b, gas,
                        sph, grw, channels_frac=ch, r_avg=r_avg))
    # multi-year: desalinated + gassy
    for z in ZS:
        T = temperature(z, T_top=-20.0)
        S = 1.5                                   # drained, near-fresh
        phi_b = phi_brine(S, T)
        gas = 0.07 - 0.04 * z                     # gas-rich, more near surface (0.07->0.03)
        sph, grw = 0.88 - 0.20*z, 0.35 + 0.25*z   # rounder, less oriented than FY
        r_avg = 0.030 + 0.012*z
        ch = 0.0                                   # too little brine to percolate
        rows.append(row(f'FYMY_my_z{int(z*100):02d}', E_matrix(T), phi_b, gas,
                        sph, grw, channels_frac=ch, r_avg=r_avg))
    write('rve_fymy.csv', rows)

# ---------------------------------------------------------------------------
# Study 7: FULL 6x6 TENSOR of the channelled base. Realistic FY C-shape base
# slices (z65..z95), full_tensor=Yes so shear anisotropy G_xz/G_yz vs G_xy is
# resolved (complete transverse-isotropy tensor, not just E_z/E_x).
def study_basetensor():
    S_cshape = {0.65: 4.5, 0.75: 5.0, 0.85: 6.0, 0.95: 8.0}
    rows = []
    for z, S in S_cshape.items():
        T = temperature(z, T_top=-20.0)
        phi_b = phi_brine(S, T)
        gas = 0.012 + 0.008*(1-z)
        sph, grw = 0.86 - 0.26*z, 0.40 + 0.32*z
        r_avg = 0.030 + 0.016*z
        ch = 0.40 if phi_b > 0.05 else 0.0
        r = row(f'BTEN_z{int(z*100):02d}', E_matrix(T), phi_b, gas, sph, grw,
                channels_frac=ch, r_avg=r_avg)
        r['full_tensor'] = 'Yes'                  # 6 load cases
        rows.append(r)
    write('rve_basetensor.csv', rows)

# ---------------------------------------------------------------------------
# Study 9 (2026-07-05): FULL 6x6 TENSOR down the WHOLE column, not just the
# base 4 slices. Same FY C-shape physics as study_basetensor(), extended to all
# ten depths z05..z95 with full_tensor=Yes, so the complete depth evolution of
# E_x/E_y/E_z, G_ij and the transverse-isotropy constants is resolved. The four
# base slices (z65..z95) reproduce the earlier basetensor microstructure exactly
# (identical formulas), so the campaign is a self-contained superset.
def study_coltensor():
    S_cshape = [7.0, 5.5, 4.8, 4.5, 4.3, 4.3, 4.5, 5.0, 6.0, 8.0]  # Cox & Weeks 1983
    rows = []
    for z, S in zip(ZS, S_cshape):
        T = temperature(z, T_top=-20.0)
        phi_b = phi_brine(S, T)
        gas = 0.012 + 0.008*(1-z)
        sph, grw = 0.86 - 0.26*z, 0.40 + 0.32*z
        r_avg = 0.030 + 0.016*z
        ch = 0.40 if phi_b > 0.05 else 0.0
        r = row(f'CTEN_z{int(z*100):02d}', E_matrix(T), phi_b, gas, sph, grw,
                channels_frac=ch, r_avg=r_avg)
        r['full_tensor'] = 'Yes'                  # 6 load cases
        rows.append(r)
    write('rve_coltensor.csv', rows)

# ---------------------------------------------------------------------------
# Study 8: STEEPER monotonic salinity to close the residual Marchenko gap.
# The first monotonic run was 15-26% too stiff in the interior because its
# salinity gradient was too gentle; here S rises faster (concave) so phi_b and
# the E knockdown grow more quickly down-column.
def study_mono2():
    rows = []
    for z in ZS:
        T = temperature(z, T_top=-20.0)
        S = 4.0 + 9.0 * z**0.6            # steeper/concave: 4.7 top -> 12.8 base
        phi_b = min(phi_brine(S, T), 0.16)  # cap: keep base physical & mesh-safe
        gas = 0.012
        sph = 0.86 - 0.26*z
        grw = 0.40 + 0.32*z
        r_avg = 0.030 + 0.016*z
        ch = 0.45 if phi_b > 0.05 else 0.0
        rows.append(row(f'MON2_z{int(z*100):02d}', E_matrix(T, scale=BEAM_FACTOR),
                        phi_b, gas, sph, grw, channels_frac=ch, r_avg=r_avg))
    write('rve_mono2.csv', rows)

if __name__ == '__main__':
    study_channel()
    study_fymy()
    study_basetensor()
    study_coltensor()
    study_mono2()
